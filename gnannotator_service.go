package genome_nexus_annotator_go

import (
	"context"
	"fmt"
	"strconv"

	gnapi "github.com/genome-nexus/genome-nexus-go-api-client/genome-nexus-public-api"
	tt "github.mskcc.org/cdsi/cdsi-protobuf/tempo/generated/v2/go"
)

type GNAnnotator interface {
	GetGenomeNexusInfo() (*gnapi.AggregateSourceInfo, error)
	AnnotateTempoMessageEvents(isoformOverrideSource string, tm *tt.TempoMessage) error
}

type GNAnnotatorService struct {
	client         *gnapi.APIClient
	ctxAccessToken context.Context
	token          string
}

func NewGNAnnotatorService(ctx context.Context, token, gnURL string) (GNAnnotator, error) {
	if len(gnURL) == 0 {
		return nil, fmt.Errorf("gnURL: %q needs to be valid", gnURL)
	}
	cfg := gnapi.NewConfiguration()
	cfg.Servers = gnapi.ServerConfigurations{
		{
			URL:         gnURL,
			Description: "Genome Nexus Annotation Server",
		},
	}
	client := gnapi.NewAPIClient(cfg)
	if len(token) > 0 {
		ctx = context.WithValue(ctx, gnapi.ContextAccessToken, token)
	}
	return GNAnnotatorService{client: client, ctxAccessToken: ctx, token: token}, nil

}

func (gn GNAnnotatorService) GetGenomeNexusInfo() (*gnapi.AggregateSourceInfo, error) {
	resp, _, err := gn.client.InfoControllerApi.FetchVersionGET(gn.ctxAccessToken).Execute()
	if err != nil {
		return resp, fmt.Errorf("Genome Nexus request failed: %v", err)
	}
	return resp, nil
}

func (gn GNAnnotatorService) AnnotateTempoMessageEvents(
	isoformOverrideSource string,
	tm *tt.TempoMessage,
) error {
	// prepare genomic locations for annotation
	genomicLocations := make([]gnapi.GenomicLocation, 0)
    for _, a := range tm.Events {
		loc := gn.getGenomicLocation(a)
		genomicLocations = append(genomicLocations, loc)
	}

	variantAnnotations, err := gn.getVariantAnnotations(isoformOverrideSource, genomicLocations)
	if err != nil {
		return err
	}

    // Build a mapping from genomic location key -> indices of records
    genomicLocationToRecordIndices := make(map[string][]int)
    for i, gl := range genomicLocations {
        key := buildGenomicLocationKey(gl)
        genomicLocationToRecordIndices[key] = append(genomicLocationToRecordIndices[key], i)
    }

    // Track which records were annotated by the response
    annotated := make([]bool, len(genomicLocations))

    // Map each returned variant annotation to the correct record(s) by key
    for _, variantAnnotation := range variantAnnotations {
        // Prefer the original variant query key when available
        var key string
        if variantAnnotation.OriginalVariantQuery != nil && *variantAnnotation.OriginalVariantQuery != "" {
            key = *variantAnnotation.OriginalVariantQuery
        } else if variantAnnotation.AnnotationSummary != nil && variantAnnotation.AnnotationSummary.GenomicLocation != nil {
            key = buildGenomicLocationKey(*variantAnnotation.AnnotationSummary.GenomicLocation)
        } else {
            continue
        }
        if indices, ok := genomicLocationToRecordIndices[key]; ok {
            for _, idx := range indices {
                gn.mapResponseToEvent(variantAnnotation, genomicLocations[idx], tm.Events[idx])
                annotated[idx] = true
            }
        }
    }

    // Any records not annotated by Genome Nexus response should be marked as failure
    for i := range genomicLocations {
        if !annotated[i] {
            tm.Events[i].AnnotationStatus = fmt.Sprintf(
                "FAILURE: No variant annotation returned for genomicLocation: %s",
                buildGenomicLocationKey(genomicLocations[i]),
            )
        }
    }
	return nil
}

func (gn GNAnnotatorService) getGenomicLocation(e *tt.Event) gnapi.GenomicLocation {
	start, _ := strconv.ParseInt(e.StartPosition, 10, 32)
	end, _ := strconv.ParseInt(e.EndPosition, 10, 32)

	gloc := gnapi.NewGenomicLocation(e.Chromosome,
		int32(start),
		int32(end),
		e.ReferenceAllele,
		resolveTumorSeqAlleleFromInput(e.ReferenceAllele, e.TumorSeqAllele1, e.TumorSeqAllele2))
	return *gloc
}

func (gn GNAnnotatorService) getVariantAnnotations(
	isoformOverrideSource string,
	genomicLocations []gnapi.GenomicLocation,
) ([]gnapi.VariantAnnotation, error) {
	fields := make([]string, 0)
	fields = append(fields, "annotation_summary")
	x := gn.client.AnnotationControllerApi.FetchVariantAnnotationByGenomicLocationPOST(gn.ctxAccessToken).
		GenomicLocations(genomicLocations)
	x = x.Fields(fields)
	x = x.IsoformOverrideSource(isoformOverrideSource)
	x = x.Token(gn.token)

	variantAnnotations, r, err := x.Execute()
	if err != nil {
		return variantAnnotations, fmt.Errorf(
			"Error calling Genome Nexus annotation service: %v\nFull HTTP response: %v\n",
			err,
			r,
		)
	}
	return variantAnnotations, nil
}

func (gn GNAnnotatorService) mapResponseToEvent(
	variantAnnotation gnapi.VariantAnnotation,
	genomicLocation gnapi.GenomicLocation,
	event *tt.Event,
) {
	if !*variantAnnotation.SuccessfullyAnnotated {
		event.AnnotationStatus = fmt.Sprintf(
			"FAILURE: Unsuccessful variant annotation for genomicLocation: %v",
			variantAnnotation,
		)
		return
	}
	canonicalTranscript := getCanonicalTranscript(variantAnnotation)
	// Taken from GN response (default)
	event.Chromosome = resolveChromosome(
		variantAnnotation,
		genomicLocation,
	) // annotationUtil.resolveChromosome(gnResponse, mRecord)
	event.StartPosition = resolveStart(variantAnnotation, genomicLocation)
	event.EndPosition = resolveEnd(
		variantAnnotation,
		genomicLocation,
	) // annotationUtil.resolveEnd(gnResponse, mRecord)
	event.Strand = resolveStrandSign(
		variantAnnotation,
	) // annotationUtil.resolveStrandSign(gnResponse, mRecord)
	event.NcbiBuild = resolveAssemblyName(
		variantAnnotation,
	) // annotationUtil.resolveAssemblyName(gnResponse, mRecord)
	event.HugoSymbol = resolveHugoSymbol(canonicalTranscript)
	event.EntrezGeneId = resolveEntrezGeneId(canonicalTranscript)
	event.VariantClassification = resolveVariantClassification(
		canonicalTranscript,
		*event,
	) // annotationUtil.resolveVariantClassification(gnResponse, canonicalTranscript, mRecord)
	event.VariantType = resolveVariantType(
		variantAnnotation,
	) // annotationUtil.resolveVariantType(gnResponse)
	event.DbsnpRs = resolveDbSnpRs(
		variantAnnotation,
	) // 	annotationUtil.resolveDbSnpRs(gnResponse, mRecord)
	event.Hgvsc = resolveHgvsc(
		canonicalTranscript,
	) // annotationUtil.resolveHgvsc(canonicalTranscript)
	event.Hgvsp = resolveHgvsp(
		canonicalTranscript,
	) // annotationUtil.resolveHgvsp(canonicalTranscript)
	event.HgvspShort = resolveHgvspShort(
		canonicalTranscript,
	) // annotationUtil.resolveHgvspShort(canonicalTranscript)
	event.TranscriptId = resolveTranscriptId(
		canonicalTranscript,
	) // annotationUtil.resolveTranscriptId(canonicalTranscript)
	event.Refseq = resolveRefSeq(
		canonicalTranscript,
	) // annotationUtil.resolveRefSeq(canonicalTranscript)
	event.Codons = resolveCodonChange(canonicalTranscript)
	event.Consequence = resolveConsequence(canonicalTranscript)
	event.ProteinPosition = resolveProteinPosition(canonicalTranscript)
	event.ExonNumber = resolveExon(canonicalTranscript)

	// Taken from GN response (Polyphen)
	event.PolyphenPrediction = resolvePolyphenPrediction(
		canonicalTranscript,
	) // annotationUtil.resolvePolyphenPrediction(canonicalTranscript)
	event.PolyphenScore = resolvePolyphenScore(
		canonicalTranscript,
	) // annotationUtil.resolvePolyphenScore(canonicalTranscript)

	// Taken from GN response (SIFT)
	event.SiftPrediction = resolveSiftPrediction(
		canonicalTranscript,
	) // annotationUtil.resolveSiftPrediction(canonicalTranscript)
	event.SiftScore = resolveSiftScore(
		canonicalTranscript,
	) // annotationUtil.resolveSiftScore(canonicalTranscript)

	// ======================================
	// separate special function needed
	event.ReferenceAllele, event.TumorSeqAllele1, event.TumorSeqAllele2 = resolveRefAndTumorSeqAlleles(
		variantAnnotation,
		*event,
		"true",
	)
	// ======================================
	event.AnnotationStatus = "SUCCESS"
}

// buildGenomicLocationKey constructs a stable key for a genomic location using
// chromosome, start, end, referenceAllele and variantAllele. This is used to
// correlate Genome Nexus responses back to the original records regardless of order.
func buildGenomicLocationKey(gl gnapi.GenomicLocation) string {
    // Chromosome
    chromosome := gl.Chromosome
    if c, ok := gl.GetChromosomeOk(); ok && c != nil && *c != "" {
        chromosome = *c
    }

    // Start
    start := gl.Start
    if s, ok := gl.GetStartOk(); ok && s != nil {
        start = *s
    }

    // End
    end := gl.End
    if e, ok := gl.GetEndOk(); ok && e != nil {
        end = *e
    }

    // Reference allele
    referenceAllele := gl.ReferenceAllele
    if r, ok := gl.GetReferenceAlleleOk(); ok && r != nil {
        referenceAllele = *r
    }

    // Variant allele
    variantAllele := gl.VariantAllele
    if v, ok := gl.GetVariantAlleleOk(); ok && v != nil {
        variantAllele = *v
    }

    return fmt.Sprintf("%s:%d:%d:%s:%s", chromosome, start, end, referenceAllele, variantAllele)
}
