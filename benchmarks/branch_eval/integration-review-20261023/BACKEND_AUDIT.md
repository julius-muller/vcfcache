# Backend-aware review

142 archived timed rows retained; 119 primary matrix cells; 3 primary heavy cells; 20 superseded cells.

Requested arm and observed backend are distinct. Historical arm `auto` means adaptive lookup, not the new availability-based CLI engine `auto`.

| Requested arm | Effective backend | Primary cells |
| --- | --- | ---: |
| auto | adaptive | 21 |
| direct | uncached | 20 |
| legacy | legacy | 20 |
| main-legacy | legacy | 21 |
| native | adaptive | 6 |
| native | native | 15 |
| original-fused | original_fused | 19 |

## Native requests that fell back

| Recipe | Scale | Job | Observed backend | Reason |
| --- | --- | --- | --- | --- |
| fastvep | wgs | 1316 | adaptive | native memory buffer budget exceeded (indexed record-count lower bound) |
| everything_delay_00500us | KPGP-00319.GRCh38.small_variants.vcf | 1320 | adaptive | native memory buffer budget exceeded (indexed record-count lower bound) |
| vep | wes | 1315 | adaptive | UNSUPPORTED: native memory buffer budget exceeded |
| vep | wgs | 1317 | adaptive | native memory buffer budget exceeded (indexed record-count lower bound) |
| veprs | wgs | 1302 | adaptive | native memory buffer budget exceeded (indexed record-count lower bound) |
| veprs-clinical | wgs | 1324 | adaptive | native memory buffer budget exceeded (indexed record-count lower bound) |

Audits are unchanged. A baseline self-check is not independent validation. Clinical mismatches remain failures, not ignored fields. Small timing differences and slopes spanning fallback backends do not identify an engine advantage.
