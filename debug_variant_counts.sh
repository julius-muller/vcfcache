#!/bin/bash
# Debug script to check variant counts at each step

CACHED_DIR="$1"
UNCACHED_DIR="$2"

if [ -z "$CACHED_DIR" ] || [ -z "$UNCACHED_DIR" ]; then
    echo "Usage: $0 <cached_dir> <uncached_dir>"
    echo "Example: $0 ft1/cached ft1/uncached"
    exit 1
fi

echo "=========================================="
echo "Variant Count Analysis"
echo "=========================================="
echo ""

echo "Input file (should be same for both):"
echo "  Original: $(bcftools index -n ft1.bcf 2>/dev/null || echo 'N/A')"
echo ""

echo "Uncached workflow:"
if [ -f "$UNCACHED_DIR/work/annotate-nocache/"*"_normalized.bcf" ]; then
    NORM=$(ls "$UNCACHED_DIR/work/annotate-nocache/"*"_normalized.bcf" 2>/dev/null | head -1)
    echo "  After norm:  $(bcftools index -n "$NORM" 2>/dev/null || echo 'N/A')"
fi
echo "  Final output: $(bcftools index -n "$UNCACHED_DIR/"*"_vc.bcf" 2>/dev/null || echo 'N/A')"
echo ""

echo "Cached workflow:"
if [ -f "$CACHED_DIR/work/annotate/"*"_normalized.bcf" ]; then
    NORM=$(ls "$CACHED_DIR/work/annotate/"*"_normalized.bcf" 2>/dev/null | head -1)
    echo "  After norm:  $(bcftools index -n "$NORM" 2>/dev/null || echo 'N/A')"
fi
if [ -f "$CACHED_DIR/work/annotate/"*"_isec_vc.bcf" ]; then
    ISEC=$(ls "$CACHED_DIR/work/annotate/"*"_isec_vc.bcf" 2>/dev/null | head -1)
    echo "  Step 1 (with cache): $(bcftools index -n "$ISEC" 2>/dev/null || echo 'N/A')"
fi
if [ -f "$CACHED_DIR/work/annotate/"*"_isec_vc_miss.bcf" ]; then
    MISS=$(ls "$CACHED_DIR/work/annotate/"*"_isec_vc_miss.bcf" 2>/dev/null | head -1)
    echo "  Step 2 (missing):   $(bcftools index -n "$MISS" 2>/dev/null || echo 'N/A')"
fi
if [ -f "$CACHED_DIR/work/annotate/"*"_missing_annotated.bcf" ]; then
    ANNO=$(ls "$CACHED_DIR/work/annotate/"*"_missing_annotated.bcf" 2>/dev/null | head -1)
    echo "  Step 3 (annotated): $(bcftools index -n "$ANNO" 2>/dev/null || echo 'N/A')"
fi
echo "  Final output: $(bcftools index -n "$CACHED_DIR/"*"_vc.bcf" 2>/dev/null || echo 'N/A')"
echo ""

echo "Difference: $(($(bcftools index -n "$CACHED_DIR/"*"_vc.bcf" 2>/dev/null || echo 0) - $(bcftools index -n "$UNCACHED_DIR/"*"_vc.bcf" 2>/dev/null || echo 0))) variants"
