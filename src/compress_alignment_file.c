#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <argp.h>
#include <htslib/sam.h>
#include <sys/types.h>
#include "data_structure_definitions.h"
#include "function_definitions.h"

/*
 * Set up the argument parser
 */
const char *argp_program_version = "abridge compress_alignment_file 1.0.0";
const char *argp_program_bug_address = "sagnikbanerjee15@gmail.com";
static char doc[] = "compress_alignment_file will accept an alignment file in SAM/BAM format and remove all redundant information. It will preserve only the information that has been requested by the user.";
static char args_doc[] = ""; // No standard arguments
							 // (i.e. arguments without "names")

static struct argp_option options[] =
{
/*
 * Options.  Field 1 in ARGP.
 * Order of fields: {NAME, KEY, ARG, FLAGS, DOC, GROUP}.
 */
{ "input_alignment_filename" , 'i' , "ALIGNMENT_FILENAME" , 0 , "Enter the name of the alignment file to be compressed" , 0 } ,
{ "input_alignment_file_format" , 'j' , "ALIGNMENT_FILE_FORMAT" , 0 , "Enter the format of the alignment file. Must be either SAM or BAM" , 0 } ,
{ "output_abridge_filename" , 'o' , "TEXT_FILENAME" , 0 , "Enter the name of the compressed file (please note that this is not the final compressed file)" , 0 } ,
{ "reference_filename" , 'g' , "reference_filename" , 0 , "Enter the name of the genome file in fasta format" , 0 } ,
{ "unmapped_filename" , 'u' , "UNMAPPED_READS_FILENAME" , 0 , "Enter the name of the file where the unmapped reads will be stored" , 0 } ,
{ "max_read_length" , 'c' , "MAX_READ_LENGTH" , 0 , "Maximum read length of short reads" , 0 } ,
{ "name_of_file_with_quality_scores" , 'q' , "QUALITY_SCORES_FILENAME" , 0 , "Enter the name of the file where the quality scores will be stored. This file will be compressed later" , 0 } ,
{ "name_of_file_with_read_names_to_short_read_names_and_NH" , 'r' , "SHORT_NAMES_NH_FILENAME" , 0 , "Enter the name of the file that contains the mapping between the long name to the short name and the NH values" , 0 } ,
{ "ended" , 'k' , "PE_OR_SE" , 0 , "Enter SE or PE to indicate whether the alignment file is SE or PE" , 0 } ,
{ "AS_tag_presence" , 'l' , "AS_tag_presence" , 0 , "Enter 1 or 0 depending on whether the AS tag is present or not" , 0 } ,

{ "flag_ignore_soft_clippings" , 's' , 0 , 0 , "Set this flag to ignore soft clippings" , 0 } ,
{ "flag_ignore_mismatches" , 'm' , 0 , 0 , "Set this flag to ignore mismatches" , 0 } ,
{ "flag_ignore_all_quality_scores" , 'p' , 0 , 0 , "Set this flag to ignore quality scores for mismatched bases and soft clips" , 0 } ,
{ "flag_ignore_unmapped_sequences" , 'e' , 0 , 0 , "Set this flag to ignore unmapped sequences along with their quality scores" , 0 } ,
{ "flag_ignore_quality_scores_for_matched_bases" , 'b' , 0 , 0 , "Set this flag to ignore quality scores for nucleotide bases that match to the provided reference" , 0 } ,
{ "flag_ignore_alignment_scores" , 'a' , 0 , 0 , "Set this flag to ignore the alignment scores (Column 5 of SAM file)" , 0 } ,
{ "skip_shortening_read_names" , 'f' , 0 , 0 , "Set this flag to skip shortening read names" , 0 } ,
{ "run_diagnostics" , 'd' , 0 , 0 , "Set this flag to run diagnostics and print out a verbose report" , 0 } ,

{ "max_reads_in_a_single_nucl_loc" , 'n' , "MAX_READS_IN_ONE_NUCL" , 0 , "Enter the value of the maximum number of input reads mapped to a single nucleotide" , 0 } ,
{ 0 , 0 , 0 , 0 , 0 , 0 } // Last entry should be all zeros in all fields
};

struct arguments
{
	/* Used by main to communicate with parse_opt. */
	// char *args[0];   // No standard arguments (without flags)
	char *input_alignment_filename;
	char *input_alignment_file_format;
	char *output_abridge_filename;
	char *reference_filename;
	char *ended;
	char *unmapped_filename;
	char *name_of_file_with_quality_scores;
	char *name_of_file_with_read_names_to_short_read_names_and_NH;
	unsigned short int flag_ignore_soft_clippings;
	unsigned short int flag_ignore_mismatches;
	unsigned short int flag_ignore_all_quality_scores;
	unsigned short int flag_ignore_unmapped_sequences;
	unsigned short int run_diagnostics;
	unsigned short int flag_ignore_quality_scores_for_matched_bases;
	unsigned short int flag_ignore_alignment_scores;
	unsigned short int skip_shortening_read_names;
	unsigned short int AS_tag_presence;
	unsigned int max_read_length;
	unsigned long long int max_reads_in_a_single_nucl_loc;
};

static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
	/*
	 * Parser. Field 2 in ARGP.
	 * Order of parameters: KEY, ARG, STATE.
	 * Parse a single option.
	 */
	/* Get the input argument from argp_parse, which we
	 know is a pointer to our arguments structure. */
	struct arguments *arguments = state->input;
	char *eptr;

	// Figure out which option we are parsing, and decide how to store it
	switch ( key )
	{
		case 'a':
			arguments->flag_ignore_alignment_scores = 1;
			break;
		case 'b':
			arguments->flag_ignore_quality_scores_for_matched_bases = 1;
			break;
		case 'c':
			arguments->max_read_length = convertStringToUnsignedInteger (arg);
			break;
		case 'd':
			arguments->run_diagnostics = 1;
			break;
		case 'e':
			arguments->flag_ignore_unmapped_sequences = 1;
			break;
		case 'f':
			arguments->skip_shortening_read_names = 1;
			break;
		case 'g':
			arguments->reference_filename = arg;
			break;
		case 'i':
			arguments->input_alignment_filename = arg;
			break;
		case 'j':
			arguments->input_alignment_file_format = arg;
			break;
		case 'k':
			arguments->ended = arg;
			break;
		case 'l':
			arguments->AS_tag_presence = convertStringToUnsignedInteger (arg);
			break;
		case 'm':
			arguments->flag_ignore_mismatches = 1;
			break;
		case 'n':
			arguments->max_reads_in_a_single_nucl_loc = convertStringToUnsignedInteger (arg);
			break;
		case 'o':
			arguments->output_abridge_filename = arg;
			break;
		case 'p':
			arguments->flag_ignore_all_quality_scores = 1;
			break;
		case 'q':
			arguments->name_of_file_with_quality_scores = arg;
			break;
		case 'r':
			arguments->name_of_file_with_read_names_to_short_read_names_and_NH = arg;
			break;
		case 's':
			arguments->flag_ignore_soft_clippings = 1;
			break;
		case 'u':
			arguments->unmapped_filename = arg;
			break;

		case ARGP_KEY_END:
			// Reached the last key.
			// Check if our REQUIRED "options" have been set to non-default values
			if ( strcmp (arguments->input_alignment_filename , "") == 0 || strcmp (arguments->input_alignment_file_format ,
					"") == 0 || strcmp (arguments->output_abridge_filename , "") == 0 || strcmp (arguments->reference_filename ,
					"") == 0 || strcmp (arguments->unmapped_filename , "") == 0 || strcmp (arguments->name_of_file_with_quality_scores ,
					"") == 0 || strcmp (arguments->name_of_file_with_read_names_to_short_read_names_and_NH ,
					"") == 0 || arguments->max_reads_in_a_single_nucl_loc == 0 )
			{
				argp_usage (state);
			}
			break;

		default:
			return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

// Our argp parser.
static struct argp argp =
{ options , parse_opt , args_doc , doc , 0 , 0 , 0 };

void compressAlignmentFile (
		char *reference_filename,
		char *input_alignment_filename,
		char *input_alignment_file_format,
		char *output_abridgefilename,
		char *ended,
		char *unmapped_filename,
		char *name_of_file_with_quality_scores,
		char *name_of_file_with_read_names_to_short_read_names_and_NH,

		unsigned short int flag_ignore_alignment_scores,
		unsigned short int flag_ignore_soft_clippings,
		unsigned short int flag_ignore_mismatches,
		unsigned short int flag_ignore_all_quality_scores,
		unsigned short int flag_ignore_unmapped_sequences,
		unsigned short int flag_ignore_quality_scores_for_matched_bases,
		unsigned short int run_diagnostics,
		unsigned short int max_reads_in_a_single_nucl_loc,
		unsigned short int skip_shortening_read_names,
		unsigned short int AS_tag_presence,
		unsigned int max_read_length)
{
	/****************************************************************************************************************************************
	 * Variable declaration
	 ****************************************************************************************************************************************/
	FILE *fhr;
	FILE *fhw_compressed;
	FILE *fhw_unmapped;
	FILE *fhw_qual;
	FILE *fhr_name_of_file_with_read_names_to_short_read_names_and_NH;

	char **split_on_tab; // List of strings to store each element of a single alignment
	char **split_on_colon; // List of strings to store tag information
	char *temp; //Useless
	char *line = NULL; // for reading each line
	char *line_name_of_file_with_read_names_to_short_read_names_and_NH = NULL;
	char *entry_in_output_file; //entry in output file
	char *prev_reference_name;
	char *curr_reference_name;
	char *reference_id_quick_read;
	char *samflag_quick_read;
	char **modified_icigars;
	char *line_to_be_written_to_file;
	char *list_of_read_names;
	char *list_of_qual_scores;
	char *qual_for_writeToFile;
	char str[100];

	size_t len = 0;
	ssize_t line_len;

	short *already_processed;

	int flag;
	int i, j, k; // Required in loops
	int number_of_tags;
	int number_of_fields; // Number of fields in each sam alignment entry
	int sam_tag_index;
	int tab_number;
	int num_items_in_alignment_pool = 0; // Items in pool
	int samflag_quick_read_index = 0;
	int compressed_ds_pool_index = 0;
	int quality_score_index = 0;
	int number_of_reference_sequences = 0;
	int reference_sequence_index = 0;

	long long int relative_position_to_previous_read_cluster;
	long long int previous_position = -1;
	long long int current_position;
	long long int number_of_records_written = 0;
	long long int number_of_records_read = 0;
	long long int num_pools_written = 0;
	long long int max_commas = 0;
	long long int curr_commas = 0;
	unsigned long long int sam_alignment_instance_pool_size = 0;
	unsigned long long int sam_alignment_instance_pool_index = 0;

	/* Variables if BAM file is provided*/
	samFile *fp_in;            // File pointer if BAM file provided
	bam_hdr_t *bamHdr;         // read header
	bam1_t *aln = bam_init1 (); // initialize an alignment

	struct Sam_Alignment *previous_alignment;
	struct Sam_Alignment *current_alignment;
	struct Sam_Alignment *sam_alignment_instance_diagnostics;
	struct Sam_Alignment *temp_alignment;
	struct Sam_Alignment **alignment_pool_same_position;
	struct Sam_Alignment **sam_alignment_instance_pool;
	struct Reference_Sequence_Info **reference_info;
	struct Whole_Genome_Sequence *whole_genome;
	struct Cigar_Items **cigar_items_instance;
	/****************************************************************************************************************************************/

	/****************************************************************************************************************************************
	 * Variable initialization
	 ****************************************************************************************************************************************/
	if ( strcmp (input_alignment_file_format , "SAM") == 0 )
	{
		fhr = fopen (input_alignment_filename , "r");
		if ( fhr == NULL )
		{
			printf ("Error! File %s not found" , input_alignment_filename);
			exit (1);
		}
	}
	else if ( strcmp (input_alignment_file_format , "BAM") == 0 )
	{
		fp_in = hts_open (input_alignment_filename , "r");
		bamHdr = sam_hdr_read (fp_in);
	}

	fhw_compressed = fopen (output_abridgefilename , "w");
	if ( fhw_compressed == NULL )
	{
		printf ("%s File cannot be created" , output_abridgefilename);
		exit (1);
	}
	fhw_unmapped = fopen (unmapped_filename , "w");
	if ( fhw_unmapped == NULL )
	{
		printf ("%s File cannot be created" , unmapped_filename);
		exit (1);
	}

	fhw_qual = fopen (name_of_file_with_quality_scores , "w");
	if ( fhw_qual == NULL )
	{
		printf ("%s File cannot be created" , name_of_file_with_quality_scores);
		exit (1);
	}
	fhr_name_of_file_with_read_names_to_short_read_names_and_NH = fopen (name_of_file_with_read_names_to_short_read_names_and_NH ,
			"r");
	if ( fhr_name_of_file_with_read_names_to_short_read_names_and_NH == NULL )
	{
		printf ("Error! File %s not found" ,
				name_of_file_with_read_names_to_short_read_names_and_NH);
		exit (1);
	}

	split_on_tab = ( char** ) malloc (sizeof(char*) * ONE_HUNDRED);
	for ( i = 0 ; i < ONE_HUNDRED ; i++ )
		split_on_tab[i] = ( char* ) malloc (sizeof(char) * ( max_read_length + TEN ));

	split_on_colon = ( char** ) malloc (sizeof(char*) * TEN);
	for ( i = 0 ; i < TEN ; i++ )
		split_on_colon[i] = ( char* ) malloc (sizeof(char) * ONE_HUNDRED);

	already_processed = ( short* ) malloc (sizeof(short) * max_reads_in_a_single_nucl_loc);
	max_reads_in_a_single_nucl_loc += 5;

	line_to_be_written_to_file = ( char* ) malloc (sizeof(char) * MAX_LINE_TO_BE_WRITTEN_TO_FILE);
	qual_for_writeToFile = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);

	reference_id_quick_read = ( char* ) malloc (sizeof(char) * ONE_THOUSAND);
	samflag_quick_read = ( char* ) malloc (sizeof(char) * ONE_THOUSAND);
	prev_reference_name = ( char* ) malloc (sizeof(char) * ONE_THOUSAND);
	prev_reference_name[0] = '\0';
	curr_reference_name = ( char* ) malloc (sizeof(char) * ONE_THOUSAND);
	curr_reference_name[0] = '\0';

	current_alignment = allocateMemorySam_Alignment (max_read_length);
	previous_alignment = allocateMemorySam_Alignment (max_read_length);
	temp_alignment = allocateMemorySam_Alignment (max_read_length);
	cigar_items_instance = ( struct Cigar_Items** ) malloc (sizeof(struct Cigar_Items*) * max_read_length);
	for ( i = 0 ; i < max_read_length ; i++ )
		cigar_items_instance[i] = allocateMemoryCigar_Items ();
	sam_alignment_instance_diagnostics = allocateMemorySam_Alignment (max_read_length);
	reference_info = ( struct Reference_Sequence_Info** ) malloc (sizeof(struct Reference_Sequence_Info*) * MAX_REFERENCE_SEQUENCES);
	for ( i = 0 ; i < MAX_REFERENCE_SEQUENCES ; i++ )
		reference_info[i] = allocateMemoryReference_Sequence_Info ();

	modified_icigars = ( char** ) malloc (sizeof(char*) * max_reads_in_a_single_nucl_loc);
	for ( i = 0 ; i < max_reads_in_a_single_nucl_loc ; i++ )
		modified_icigars[i] = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);

	sam_alignment_instance_pool = ( struct Sam_Alignment* ) malloc (sizeof(struct Sam_Alignment*) * max_reads_in_a_single_nucl_loc);
	for ( i = 0 ; i < max_reads_in_a_single_nucl_loc ; i++ )
		sam_alignment_instance_pool[i] = allocateMemorySam_Alignment (max_read_length);
	sam_alignment_instance_pool_size = max_reads_in_a_single_nucl_loc;
	/****************************************************************************************************************************************/

	/*
	 * Write the first line in output file
	 */
	line_to_be_written_to_file[0] = '\0';
	strcat(line_to_be_written_to_file , "flag_ignore_mismatches:");
	convertUnsignedIntegerToString (str , flag_ignore_mismatches);
	strcat(line_to_be_written_to_file , str);
	strcat(line_to_be_written_to_file , "\t");
	strcat(line_to_be_written_to_file , "flag_ignore_soft_clippings:");
	convertUnsignedIntegerToString (str , flag_ignore_soft_clippings);
	strcat(line_to_be_written_to_file , str);
	strcat(line_to_be_written_to_file , "\t");
	strcat(line_to_be_written_to_file , "flag_ignore_unmapped_sequences:");
	convertUnsignedIntegerToString (str , flag_ignore_unmapped_sequences);
	strcat(line_to_be_written_to_file , str);
	strcat(line_to_be_written_to_file , "\t");
	strcat(line_to_be_written_to_file , "flag_ignore_all_quality_scores:");
	convertUnsignedIntegerToString (str , flag_ignore_all_quality_scores);
	strcat(line_to_be_written_to_file , str);
	strcat(line_to_be_written_to_file , "\t");
	strcat(line_to_be_written_to_file ,
			"flag_ignore_quality_scores_for_matched_bases:");
	convertUnsignedIntegerToString (str ,
			flag_ignore_quality_scores_for_matched_bases);
	strcat(line_to_be_written_to_file , str);
	strcat(line_to_be_written_to_file , "\t");
	strcat(line_to_be_written_to_file , "flag_ignore_alignment_scores:");
	convertUnsignedIntegerToString (str , flag_ignore_alignment_scores);
	strcat(line_to_be_written_to_file , str);
	strcat(line_to_be_written_to_file , "\t");
	strcat(line_to_be_written_to_file , "flag_skip_shortening_read_names:");
	convertUnsignedIntegerToString (str , skip_shortening_read_names);
	strcat(line_to_be_written_to_file , str);
	strcat(line_to_be_written_to_file , "\t");

	strcat(line_to_be_written_to_file , "AS_tag_presence:");
	convertUnsignedIntegerToString (str , AS_tag_presence);
	strcat(line_to_be_written_to_file , str);
	strcat(line_to_be_written_to_file , "\t");

	strcat(line_to_be_written_to_file , "max_read_length:");
	convertUnsignedIntegerToString (str , max_read_length);
	strcat(line_to_be_written_to_file , str);
	strcat(line_to_be_written_to_file , "\t");

	strcat(line_to_be_written_to_file , "max_reads_in_a_single_nucl_loc:");
	convertUnsignedIntegerToString (str , max_reads_in_a_single_nucl_loc);
	strcat(line_to_be_written_to_file , str);
	strcat(line_to_be_written_to_file , "\t");

	strcat(line_to_be_written_to_file , "\n");
	fprintf (fhw_compressed , "%s" , line_to_be_written_to_file);

	/*
	 * For SAM file advance the pointer to the first alignment
	 */

	if ( strcmp (input_alignment_file_format , "SAM") == 0 )
	{
		while ( ( line_len = getline ( &line , &len , fhr) ) != -1 )
			if ( line[0] != '@' ) break;
	}
	if ( strcmp (input_alignment_file_format , "BAM") == 0 )
	{
		sam_read1 (fp_in , bamHdr , aln);
	}
	//return;
	do
	{

		/***************************************************************************************
		 * Read a line from the short read names file
		 ****************************************************************************************/
		/*
		 if ( skip_shortening_read_names == 0 )
		 {
		 */
		getline ( &line_name_of_file_with_read_names_to_short_read_names_and_NH ,
				&len ,
				fhr_name_of_file_with_read_names_to_short_read_names_and_NH);

		splitByDelimiter (line_name_of_file_with_read_names_to_short_read_names_and_NH ,
				'\t' ,
				split_on_tab);
		strcpy(current_alignment->read_name , split_on_tab[3]);
		strcpy(current_alignment->NH , split_on_tab[2]);

		/*}
		 */

		/****************************************************************************************/

		/****************Collect Unmapped reads*************/
		if ( current_alignment->samflag == 4 )
		{
			if ( flag_ignore_unmapped_sequences == 0 )
			{
				//Write the unmapped reads into file
				fprintf (fhw_unmapped , "%s" , current_alignment->sequence);
				fprintf (fhw_unmapped , "%s" , "\n");
				//for ( i = 0 ; current_alignment->quality_scores[i] != '\0' ; i++ )
				// current_alignment->quality_scores[i] -= QUAL_SCORE_ADJUSTMENT;
				fprintf (fhw_qual ,
						"%s" ,
						current_alignment->quality_scores);
				fprintf (fhw_qual , "%s" , "\n");
				fprintf (fhw_qual , "%s" , "\n");
				fprintf (fhw_qual , "%s" , "\n");
			}
		}
		else
		{
			if ( strcmp (input_alignment_file_format , "SAM") == 0 )
			{
				line_len = getline ( &line , &len , fhr);
				//printf ("\nLine length %d" , line_len);
				prepareSingleRecordFromAlignmentFile (line ,
						fp_in , // File pointer if BAM file provided
						bamHdr ,		// read header
						aln ,
						fhr ,
						current_alignment ,
						ended ,
						input_alignment_file_format ,
						AS_tag_presence ,
						flag_ignore_alignment_scores ,
						flag_ignore_soft_clippings ,
						flag_ignore_mismatches ,
						flag_ignore_all_quality_scores ,
						flag_ignore_unmapped_sequences ,
						flag_ignore_quality_scores_for_matched_bases ,
						max_read_length ,
						split_on_tab ,
						split_on_colon ,
						cigar_items_instance);
			}

			else if ( strcmp (input_alignment_file_format , "BAM") == 0 )
			{
				line_len = sam_read1 (fp_in , bamHdr , aln);
				/*if ( total_number_of_alignments % 10000 == 0 )
				 printf ("\nLine length %d total_number_of_alignments %d" ,
				 line_len ,
				 total_number_of_alignments);
				 */
				//fflush (stdout);
			}
		}

		if ( line_len <= 0 ) break;
	} while ( 1 );

}

int main (int argc, char *argv[])
{
	/****************************************************************************************************************************************
	 * Named CLI
	 ****************************************************************************************************************************************/
	struct arguments arguments;

	// Parse our arguments; every option seen by parse_opt will be reflected in arguments.
	// Default values.
	arguments.input_alignment_filename = ""; // Empty string - only contains null character
	arguments.input_alignment_file_format = "";
	arguments.output_abridge_filename = "";
	arguments.reference_filename = "";
	arguments.unmapped_filename = "";
	arguments.name_of_file_with_quality_scores = "";
	arguments.ended = "";
	arguments.name_of_file_with_read_names_to_short_read_names_and_NH = "";
	arguments.flag_ignore_soft_clippings = 0;
	arguments.flag_ignore_mismatches = 0;
	arguments.flag_ignore_all_quality_scores = 0;
	arguments.flag_ignore_unmapped_sequences = 0;
	arguments.flag_ignore_quality_scores_for_matched_bases = 0;
	arguments.flag_ignore_alignment_scores = 0;
	arguments.run_diagnostics = 0;
	arguments.max_reads_in_a_single_nucl_loc = 0;
	arguments.skip_shortening_read_names = 0;

	argp_parse ( &argp , argc , argv , 0 , 0 , &arguments);
	/****************************************************************************************************************************************
	 * Variable declaration
	 ****************************************************************************************************************************************/
	char input_alignment_filename[MAX_FILENAME_LENGTH];
	char input_alignment_file_format[TEN];
	char output_abridgefilename[MAX_FILENAME_LENGTH];
	char reference_filename[MAX_FILENAME_LENGTH];
	char unmapped_filename[MAX_FILENAME_LENGTH];
	char name_of_file_with_quality_scores[MAX_FILENAME_LENGTH];
	char name_of_file_with_read_names_to_short_read_names_and_NH[MAX_FILENAME_LENGTH];
	char ended[TEN];

	unsigned short int flag_ignore_soft_clippings;
	unsigned short int flag_ignore_mismatches;
	unsigned short int flag_ignore_all_quality_scores;
	unsigned short int flag_ignore_unmapped_sequences;
	unsigned short int run_diagnostics;
	unsigned short int flag_ignore_quality_scores_for_matched_bases;
	unsigned short int flag_ignore_alignment_scores;
	unsigned short int skip_shortening_read_names;
	unsigned short int AS_tag_presence;

	unsigned int max_read_length;
	unsigned long long int max_reads_in_a_single_nucl_loc;
	/****************************************************************************************************************************************/

	/****************************************************************************************************************************************
	 * Variable initialization
	 ****************************************************************************************************************************************/
	strcpy(reference_filename , arguments.reference_filename);
	strcpy(input_alignment_filename , arguments.input_alignment_filename);
	strcpy(input_alignment_file_format , arguments.input_alignment_file_format);
	strcpy(output_abridgefilename , arguments.output_abridge_filename);
	strcpy(ended , arguments.ended);
	strcpy(unmapped_filename , arguments.unmapped_filename);
	strcpy(name_of_file_with_quality_scores ,
			arguments.name_of_file_with_quality_scores);
	strcpy(name_of_file_with_read_names_to_short_read_names_and_NH ,
			arguments.name_of_file_with_read_names_to_short_read_names_and_NH);

	flag_ignore_alignment_scores = arguments.flag_ignore_alignment_scores; // Ignore the column 5 of SAM alignment file which is often set to 255 and also the AS tag if one is provided
	flag_ignore_soft_clippings = arguments.flag_ignore_soft_clippings;
	flag_ignore_mismatches = arguments.flag_ignore_mismatches;
	flag_ignore_all_quality_scores = arguments.flag_ignore_all_quality_scores;
	flag_ignore_unmapped_sequences = arguments.flag_ignore_unmapped_sequences;
	flag_ignore_quality_scores_for_matched_bases = arguments.flag_ignore_quality_scores_for_matched_bases;
	run_diagnostics = arguments.run_diagnostics;
	max_reads_in_a_single_nucl_loc = arguments.max_reads_in_a_single_nucl_loc;
	skip_shortening_read_names = arguments.skip_shortening_read_names;
	AS_tag_presence = arguments.AS_tag_presence;
	max_read_length = arguments.max_read_length;
	/****************************************************************************************************************************************/

	compressAlignmentFile (reference_filename ,
			input_alignment_filename ,
			input_alignment_file_format ,
			output_abridgefilename ,
			ended ,
			unmapped_filename ,
			name_of_file_with_quality_scores ,
			name_of_file_with_read_names_to_short_read_names_and_NH ,
			flag_ignore_alignment_scores ,
			flag_ignore_soft_clippings ,
			flag_ignore_mismatches ,
			flag_ignore_all_quality_scores ,
			flag_ignore_unmapped_sequences ,
			flag_ignore_quality_scores_for_matched_bases ,
			run_diagnostics ,
			max_reads_in_a_single_nucl_loc ,
			skip_shortening_read_names ,
			AS_tag_presence ,
			max_read_length);
	return 0;
}
