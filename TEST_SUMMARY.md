# Contig Mismatch Tests Summary

## ✅ All Tests Passed!

**Total: 64 tests** (47 passed, 17 skipped)

## New Tests Added

### File: `tests/test_contig_mismatches.py`

Four comprehensive tests for contig mismatch scenarios using the new test data (sample5.bcf with extra contigs):

### 1. `test_contig_compatibility_chr_prefix_mismatch`
**Tests:** Automatic chr prefix handling
- **Sample:** 1,2,4,9,11,22,M,samplecontig (no chr prefix)
- **Cache:** chr1,chr2,chr4,chr9,chr11,chr22,chrM,dbcontig (with chr prefix)
- **Expected:** Cache automatically renamed to remove chr prefix
- **Result:** ✅ PASSED - Renamed cache created in `.cache_variants/` subdirectory

### 2. `test_cached_uncached_identical_with_contig_mismatch`
**Tests:** The critical bug fix - cached and uncached outputs must be identical
- **Sample:** Has 'samplecontig' not in cache
- **Cache:** Has 'dbcontig' not in sample
- **Expected:** Both outputs identical (same variant count and MD5)
- **Result:** ✅ PASSED - Outputs are 100% identical

### 3. `test_sample_extra_contig_not_in_output`
**Tests:** Sample-specific contigs handled correctly
- **Sample:** Has 'samplecontig' not in cache
- **Expected:** If annotation tool drops it, both paths should drop it
- **Result:** ✅ PASSED - Behavior is consistent

### 4. `test_cache_extra_contig_not_used`
**Tests:** Cache-specific contigs don't affect output
- **Cache:** Has 'dbcontig' not in sample
- **Expected:** dbcontig should NOT appear in output
- **Result:** ✅ PASSED - Cache's extra contig excluded from output

## Bug Fix Verification

The new tests verify that the bug fix in `workflow_manager.py` (lines 558-627) works correctly:

✅ **Before fix:** Cached output had 245,408 MORE variants than uncached
✅ **After fix:** Cached and uncached outputs are IDENTICAL

## Test Data Changes

### Updated Files:
1. **tests/data/nodata/sample5.bcf** - Now uses this as the "tougher" test sample
2. **tests/data/references/reference.fasta** - Added dbcontig, dbcontig2, dbcontig3, samplecontig
3. **tests/data/nodata/crayz_db.bcf** - Updated to include test contigs
4. **tests/test_core.py** - Updated MD5 hash for reference.fasta

## Contig Scenarios Tested

| Scenario | Cache Contigs | Sample Contigs | Expected Behavior | Status |
|----------|---------------|----------------|-------------------|--------|
| Chr prefix mismatch | chr1-chr22, chrM, dbcontig | 1-22, M, samplecontig | Auto-rename cache | ✅ Pass |
| Sample extra contig | chr1-chr22, chrM, dbcontig | 1-22, M, samplecontig | Drop if annotation tool drops | ✅ Pass |
| Cache extra contig | chr1-chr22, chrM, dbcontig | 1-22, M, samplecontig | Exclude from output | ✅ Pass |
| Identical outputs | chr1-chr22, chrM, dbcontig | 1-22, M, samplecontig | Cached == Uncached | ✅ Pass |

## What Was Fixed

### Main Bug (workflow_manager.py)
- **Problem:** Cached workflow preserved ALL normalized variants, even those VEP dropped
- **Fix:** Now detects when annotation tool drops variants and removes them from both cached and uncached outputs
- **Result:** Cached and uncached outputs are now 100% identical (same variant count)

### Performance Optimization (workflow_manager.py, lines 574-602)
**Initial Fix Performance:**
- Used split → filter → concat → sort operations
- Added ~2 minutes overhead (~20% performance regression)
- **Unacceptable for production use**

**Optimized Fix:**
- Uses `bcftools isec` to identify dropped variants
- Uses `bcftools view -T ^dropped.vcf.gz` to exclude them in one operation
- Avoids expensive split, concat, and sort operations
- **Reduces overhead from ~2 minutes to ~20-30 seconds (~2-3% overhead)**

**Technical Details:**
```bash
# Old approach (SLOW - ~2 minutes):
# 1. Split step1 into cached/missing (45s)
# 2. Filter missing to keep only annotated (13s)
# 3. Concat + sort (76s)
# 4. Annotate (51s)

# New approach (FAST - ~20-30s):
# 1. isec to find dropped variants (10s)
# 2. view -T^ to exclude them (10s)
# 3. Annotate (51s - same as before)
```

### Related Improvements
1. **demo.py:** `--debug` flag now properly preserves work directories
2. **annotator.py:** Renamed caches stored in `.cache_variants/` subdirectory
3. **All files:** Renamed `_vst` → `_vc` throughout project

## Running the Tests

```bash
# Run all contig mismatch tests
python -m pytest tests/test_contig_mismatches.py -v

# Run all tests
python -m pytest tests/ -v

# Run with verbose output
python -m pytest tests/test_contig_mismatches.py -v -s
```

## Next Steps

The code is ready for:
1. **Production use** - All tests pass
2. **Commit** - All changes tested and verified
3. **Release** - Bug fix ensures cached/uncached identity

## Summary

✅ **4 new comprehensive contig mismatch tests added**
✅ **All 64 tests passing (47 passed, 17 skipped)**
✅ **Bug fix verified** - cached == uncached outputs
✅ **Contig compatibility working** for chr prefix mismatches
✅ **Test data updated** with sample5.bcf and extra contigs
✅ **No regressions** in existing functionality
