# 
import os
import re
from os.path import join
from IPython.display import HTML
import pandas as pd

# functions
def run_qc(fastq_input_folder, quality_control_output, quality_control_summary_output, default_fastq_ext='fastq', threads=20):
    """
    Runs FastQC program on all sequencing file in the given folder and generates the summarized reports of QC.
    It is assumed that FastQC and MultiQC is installed and available from the environment.

    :param fastq_input_folder:              Input folder of short read sequences. It is assumed that the extension of the files are 'default_fastq_ext' ( i.e. .fastq, .fq, etc).
    :param quality_control_output:          Output folder for individual fastq files
    :param quality_control_summary_output:  Output folder for the summarized output (MultiQC)
    :param default_fastq_ext:               Default extension of sequencing files ( i.e. .fastq, .fq, etc).
    :param threads:                         Number of CPU cores to be used during the analysis
    :return:
    """

    raw_fastqc_cmd = 'fastqc --outdir {0} \
     -t {1} \
     {2} \
    '.format(quality_control_output, threads, join(fastq_input_folder, '*' + default_fastq_ext))
    print('Running raw fastqc: %s' % raw_fastqc_cmd)
    os.system(raw_fastqc_cmd)

    # doing multiQC
    raw_multiqc_cmd = 'multiqc -f --interactive  --outdir {0} \
    {1} \
    '.format(quality_control_summary_output, join(quality_control_output, '*'))
    print('Running raw multiqc: %s' % raw_multiqc_cmd)
    os.system(raw_multiqc_cmd)

def get_qc_html_table_pair_end(output_dir):
    """
    Generates a HTML linkable output dataframe
    :param output_dir: 
    :return: 
    """
    pd.set_option('display.max_colwidth', -1)
    seq2output = get_file_pairs(output_dir, suffix_filter_to = '_fastqc.html')
    act_columns = ['ID' , 'Forward Reads', 'Reverse Reads']
    qc_results_table = pd.DataFrame.from_dict(seq2output, orient='index').reset_index()
    qc_results_table.columns = act_columns
    qc_results_table['Forward Reads'] = qc_results_table['Forward Reads'].apply(lambda x:
                                                        '<a href="https://localhost:8888/files/{0}">{1}</a>'.format(join(output_dir, x), x))
    qc_results_table['Reverse Reads'] = qc_results_table['Reverse Reads'].apply(lambda x:
                                                        '<a href="https://localhost:8888/files/{0}">{1}</a>'.format(join(output_dir, x), x))
    return qc_results_table

def get_qc_html_table_single(output_dir):
    """
    Generates a HTML linkable output dataframe
    :param output_dir:  which directory should be checked for fastqc data
    :return:
    """
    pd.set_option('display.max_colwidth', -1)
    seq2output = get_file_pairs(output_dir, suffix_filter_to='_fastqc.html')
    act_columns = ['ID', 'sample']
    qc_results_table = pd.DataFrame.from_dict(seq2output, orient='index').reset_index()




def get_illumina_pairs(save_path, filtering_string='junk', fastq_ext='.fastq'):
    """
    Assigns the sequencing pairs to an ID (sequencing)
    :param save_path:         Input folder finding the pairs in
    :param filtering_string:  Ignore those files that contain that substring. Typically 'junk' parts should be ignored
    :return:                  Dictionary where the keys are the ids and the values are the filenames of read pairs.
    """

    raw_files = sorted(
        [f for path, dfd, filenames in os.walk(save_path) for f in filenames if os.path.splitext(f)[1] == fastq_ext
         and filtering_string not in f])

    raw_prefix = list(set([get_prefix_from_filename(filename) for filename in raw_files]))
    prefix_to_pair = {}
    for filename in raw_files:
        act_prefix = get_prefix_from_filename(filename)

        if act_prefix in prefix_to_pair:
            prefix_to_pair[act_prefix].append(join(save_path, filename))
        else:
            prefix_to_pair[act_prefix] = [join(save_path, filename)]
    return prefix_to_pair

def get_prefix_from_filename(filename):
    """
    Gets the prefix (i.e. id) from a filename of a fastq file.
    :param filename:  filename of a fastq file
    :return:
    """
    # longest prefix before R1 or R2
    # only for miseq
    pair_end_match = re.compile("(\.fastq)|(\.fq)")
    match = pair_end_match.search(filename)
    try:
        prefix = filename[0:match.start()]
    except AttributeError:
        prefix = None
    return prefix

def get_trimmomatic_cmd_default_pair_end(trimmed_output, fastq_files, trimmomatic_path, trimmomatic_parameters):
    """
    Preprocesses a pair-end sequencing data given the parameters and writes the preprocessed reads into a new folder.
    It is assumed that the read pair has the same prefix and the forward one contains the substring 'R1' while the reverse one contains 'R2'.

    :param trimmed_output:          Output folder for the preprocessed reads. Preprocessed fastq files will contain the 'trimmed' substring
    :param fastq_files:             List of input fastq files (one forward and one reverse sequencing short read)
    :param trimmomatic_path:        Absolute path to the trimmomatic program
    :param trimmomatic_parameters:  Parameters for the trimmomatic
    :return:                        A runnable command
    """
    input_filename_forward = os.path.basename(fastq_files[0])
    input_filename_reverse = os.path.basename(fastq_files[1])

    (base_filename_forward, ext_forward) = os.path.splitext(input_filename_forward)
    (base_filename_reverse, ext_reverse) = os.path.splitext(input_filename_reverse)

    output_fastq_file_forward = join(trimmed_output, base_filename_forward + '.trimmed' + ext_forward)
    output_fastq_file_reverse = join(trimmed_output, base_filename_reverse + '.trimmed' + ext_reverse)

    output_fastq_junk_file_forward = join(trimmed_output, base_filename_forward + '.junk.trimmed' + ext_forward)
    output_fastq_junk_file_reverse = join(trimmed_output, base_filename_reverse + '.junk.trimmed' + ext_reverse)

    trimmomatic_cmd = '{0} PE \
                   {1} {2}                   \
                   {3} {4}                   \
                   {5} {6}                  \
                   {7} \
    '.format(trimmomatic_path, fastq_files[0], fastq_files[1],
             output_fastq_file_forward, output_fastq_junk_file_forward,
             output_fastq_file_reverse, output_fastq_junk_file_reverse,
             trimmomatic_parameters)
    return trimmomatic_cmd

def get_trimmomatic_cmd_default_single(trimmed_output, fastq_files, trimmomatic_path, trimmomatic_parameters):
    """
    Preprocesses a pair-end sequencing data given the parameters and writes the preprocessed reads into a new folder.
    It is assumed that the read pair has the same prefix and the forward one contains the substring 'R1' while the reverse one contains 'R2'.

    :param trimmed_output:          Output folder for the preprocessed reads. Preprocessed fastq files will contain the 'trimmed' substring
    :param fastq_files:             List of input fastq files (one forward and one reverse sequencing short read)
    :param trimmomatic_path:        Absolute path to the trimmomatic program
    :param trimmomatic_parameters:  Parameters for the trimmomatic
    :return:                        A runnable command
    """
    input_filename_forward = os.path.basename(fastq_files[0])
    (base_filename_forward, ext_forward) = os.path.splitext(input_filename_forward)
    output_fastq_file_forward = join(trimmed_output, base_filename_forward + '.trimmed' + ext_forward)
    trimmomatic_cmd = '{0} SE \
                   {1}                    \
                   {2}                  \
                   {3}                  \
    '.format(trimmomatic_path, fastq_files[0],
             output_fastq_file_forward,
             trimmomatic_parameters)
    return trimmomatic_cmd

def get_trimmomatic_cmd_default(trimmed_output, fastq_files, trimmomatic_path, trimmomatic_parameters):
    """
    Preprocesses a pair-end sequencing data given the parameters and writes the preprocessed reads into a new folder.
    It is assumed that the read pair has the same prefix and the forward one contains the substring 'R1' while the reverse one contains 'R2'.

    :param trimmed_output:          Output folder for the preprocessed reads. Preprocessed fastq files will contain the 'trimmed' substring
    :param fastq_files:             List of input fastq files (one forward and one reverse sequencing short read)
    :param trimmomatic_path:        Absolute path to the trimmomatic program
    :param trimmomatic_parameters:  Parameters for the trimmomatic
    :return:                        A runnable command
    """

    if len(fastq_files) == 2:
        trimmomatic_cmd = get_trimmomatic_cmd_default_pair_end(trimmed_output, fastq_files, trimmomatic_path, trimmomatic_parameters)
    elif len(fastq_files) == 1:
        trimmomatic_cmd = get_trimmomatic_cmd_default_single(trimmed_output, fastq_files, trimmomatic_path, trimmomatic_parameters)
    else:
        raise Exception('Unexpected fastq structure. Files: {0}'.format(str(fastq_files)))

    return trimmomatic_cmd
