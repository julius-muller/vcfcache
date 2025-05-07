process CaptureToolVersions {
    publishDir "${output_dir}", mode: 'copy'

    input:
    val sample_name
    path output_dir

    output:
    path "${sample_name}_final.tool_version.log", emit: version_log

    script:
    """
	echo "${params.bcftools_cmd} version:" > ${sample_name}_final.tool_version.log
    ${params.bcftools_cmd} --version >> ${sample_name}_final.tool_version.log
    """

    stub:
    """
    touch ${sample_name}_final.tool_version.log
    """
}

process ValidateInputs {
	stageInMode 'symlink'

	input:
	path vcf
	path vcf_index  // Explicitly passing the index file
	path reference
	path reference_index
	path chr_add

	output:
	val true, emit: validated
	path "validation_results.log", emit: validation_log

	script:
	"""
	# Check file existence
	test -e "${vcf}" || { echo "VCF file not found: ${vcf}"; exit 1; }
	test -e "${vcf_index}" || { echo "VCF index file not found: ${vcf_index}"; exit 1; }
	test -e "${reference}" || { echo "Reference file not found: ${reference}"; exit 1; }
	test -e "${reference_index}" || { echo "Reference index file not found: ${reference_index}"; exit 1; }
	test -e "${chr_add}" || { echo "Chr add file not found: ${chr_add}"; exit 1; }

	# Basic header check
	${params.bcftools_cmd} view -h "${vcf}" | head -n 1 >/dev/null || { echo "Invalid VCF/BCF format"; exit 1; }

    # Use the validate_vcf_ref.sh script for all validation checks
    cp \$VCFSTASH_ROOT/tools/validate_vcf_ref.sh ./
    chmod +x validate_vcf_ref.sh

    # Run the validation script with required parameters
    ./validate_vcf_ref.sh "${vcf}" "${reference}" "${chr_add}" > validation_results.log 2>&1

    # Check the exit code of the validation script
    VALIDATION_EXIT_CODE=\$?

    if [ \$VALIDATION_EXIT_CODE -ne 0 ]; then
        echo "Validation failed with exit code \$VALIDATION_EXIT_CODE" >> validation_results.log
        cat validation_results.log
        exit \$VALIDATION_EXIT_CODE
    fi

    echo "All validation checks passed successfully" >> validation_results.log


      """
  }


workflow UTILS {
    take:
    sample_name
    output_dir
    vcf
    chr_add
    reference

    main:
    // Make sure reference file index exists
    reference_index = file("${reference}.fai")

    // Explicitly look for VCF index file - try both csi and tbi formats
    vcf_index = file("${vcf}.{csi,tbi}", checkIfExists: true)[0]

    validateInputResult = ValidateInputs(
        vcf,
        vcf_index,
        reference,
        reference_index,
        chr_add,
    )

    emit:
    validate = validateInputResult.validated
    validation_log = validateInputResult.validation_log
}
