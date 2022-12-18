#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <argp.h>
#include <htslib/sam.h>
#include "data_structure_definitions.h"
#include "function_definitions.h"

// Set up the argument parser
const char *argp_program_version = "abridge compute_information_for_better_memory_management 1.1.0";
const char *argp_program_bug_address = "sagnikbanerjee15@gmail.com";
static char doc[] = "compute_information_for_better_memory_management will accept an alignment file in SAM/BAM format calculate the highest number of reads that are mapped to any nucleotide. Additionally, the program will also output the total number of alignments and the maximum read length";
static char args_doc[] = ""; // No standard arguments
							 // (i.e. arguments without "names")

static struct argp_option options[] =
{
/*
 * Options.  Field 1 in ARGP.
 * Order of fields: {NAME, KEY, ARG, FLAGS, DOC, GROUP}.
 */
{ "input_alignment_filename" , 'i' , "SAM_FILENAME" , 0 , "Enter the name of the SAM file to be compressed" , 0 } ,
{ "input_alignment_file_format" , 'j' , "ended" , 0 , "Enter the format of the alignment file. Must be either SAM or BAM" , 0 } ,
{ "summary_information_outputfilename" , 'o' , "TEXT_FILENAME" , 0 , "Enter the name of the output file that will contain (1) the value of maximum number of reads mapped to a single nucleotide (2) Total number of alignments (3) Maximum read length" , 0 } ,
{ "ended" , 'e' , "FILE_FORMAT" , 0 , "Enter whether the sample is single ended or paired ended" , 0 } ,
{ "AS_tag_presence" , 's' , "AS_TAG_PRESENSE" , 0 , "Enter 0 or 1 depending on whether the AS tag is present or not" , 0 } ,
{ 0 , 0 , 0 , 0 , 0 , 0 } // Last entry should be all zeros in all fields
};

struct arguments
{
	/* Used by main to communicate with parse_opt. */
	// char *args[0];   // No standard arguments (without flags)
	char *input_alignment_filename; // Empty string - only contains null character
	char *input_alignment_file_format;
	char *summary_information_outputfilename;
	char *ended;
	unsigned short int AS_tag_presence;
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

	// Figure out which option we are parsing, and decide how to store it
	switch ( key )
	{
		case 'i':
			arguments->input_alignment_filename = arg;
			break;
		case 'j':
			arguments->input_alignment_file_format = arg;
			break;
		case 'o':
			arguments->summary_information_outputfilename = arg;
			break;
		case 'e':
			arguments->ended = arg;
			break;
		case 's':
			arguments->AS_tag_presence = convertStringToUnsignedInteger (arg);
			break;
		case ARGP_KEY_END:
			// Reached the last key.
			// Check if our inputsamfilename and outputfilename REQUIRED "options" have been set to non-default values
			/*if ( strcmp (arguments->input_alignment_filename , "") == 0 || strcmp (strupr (arguments->input_alignment_file_format) ,
			 "SAM") != 0 || strcmp (strupr (arguments->input_alignment_file_format) ,
			 "BAM") != 0 || strcmp (arguments->summary_information_outputfilename ,
			 "") == 0 || strcmp (arguments->ended , "") == 0 )
			 {
			 argp_usage (state);
			 }*/
			break;

		default:
			return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

// Our argp parser.
static struct argp argp =
{ options , parse_opt , args_doc , doc , 0 , 0 , 0 };

void findSummaryInformation (
		char *input_alignment_filename,
		char *ended,
		char *input_alignment_file_format,
		char *summary_information_outputfilename,
		unsigned short int AS_tag_presence)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i;
	unsigned short int max_read_length, current_read_length;
	unsigned long long int total_number_of_alignments;
	unsigned long long int current_reference_position,
			previous_reference_position;
	unsigned long long int
			maximum_number_of_reads_mapped_to_a_single_reference_nucleotide,
			number_of_reads_mapped_to_the_current_reference_nucleotide;

	FILE *fhr;
	FILE *fhw;

	/* Variables if BAM file is provided*/
	samFile *fp_in;            // File pointer if BAM file provided
	bam_hdr_t *bamHdr;         // read header
	bam1_t *aln = bam_init1 (); // initialize an alignment

	size_t len = 0;
	ssize_t line_len;

	char *str;
	char *temp;        // Useless
	char *line = NULL; // for reading each line
	char **split_line; // List of strings to store each element of a single alignment
	char **split_tags;
	char *current_reference_name, *previous_reference_name;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
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
	fhw = fopen (summary_information_outputfilename , "w");
	if ( fhw == NULL )
	{
		printf ("%s File cannot be created" ,
				summary_information_outputfilename);
		exit (1);
	}

	max_read_length = 0;
	maximum_number_of_reads_mapped_to_a_single_reference_nucleotide = 0;
	number_of_reads_mapped_to_the_current_reference_nucleotide = 0;

	split_line = ( char** ) malloc (sizeof(char*) * ONE_HUNDRED);
	for ( i = 0 ; i < ONE_HUNDRED ; i++ )
		split_line[i] = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);

	split_tags = ( char** ) malloc (sizeof(char*) * ONE_HUNDRED);
	for ( i = 0 ; i < ONE_HUNDRED ; i++ )
		split_tags[i] = ( char* ) malloc (sizeof(char) * TEN);

	current_reference_name = ( char* ) malloc (sizeof(char) * ONE_HUNDRED);
	current_reference_name[0] = '\0';
	previous_reference_name = ( char* ) malloc (sizeof(char) * ONE_HUNDRED);
	previous_reference_name[0] = '\0';
	str = ( char* ) malloc (sizeof(char) * ONE_THOUSAND);

	total_number_of_alignments = 0;
	str[0] = '\0';
	/********************************************************************/

	if ( strcmp (input_alignment_file_format , "SAM") == 0 )
	{
		while ( ( line_len = getline ( &line , &len , fhr) ) != -1 )
			if ( line[0] != '@' ) break;
	}
	if ( strcmp (input_alignment_file_format , "BAM") == 0 )
	{
		sam_read1 (fp_in , bamHdr , aln);
	}

	do
	{
		total_number_of_alignments += 1;
		if ( total_number_of_alignments % 100000 == 0 )
			printf ("\nTotal number of alignments %d" ,
					total_number_of_alignments);

		if ( strcmp (input_alignment_file_format , "BAM") == 0 )
		{
			current_reference_name = bamHdr->target_name[aln->core.tid];
			current_read_length = aln->core.l_qseq;
			current_reference_position = aln->core.pos + 1;
		}
		else if ( strcmp (input_alignment_file_format , "SAM") == 0 )
		{
			splitByDelimiter (line , '\t' , split_line);
			current_read_length = strlen (split_line[9]);
			strcpy(current_reference_name , split_line[2]);
			current_reference_position = convertStringToUnsignedInteger (split_line[3]);
		}
		if ( current_reference_position == 0 ) break; //Unaligned read

		if ( max_read_length < current_read_length )
			max_read_length = current_read_length;

		if ( strlen (previous_reference_name) == 0 )
		{
			strcpy(previous_reference_name , current_reference_name);
		}
		else if ( strcmp (current_reference_name , previous_reference_name) != 0 ) // Different reference sequence encountered
		{
			previous_reference_position = -1;
		}

		//if ( current_reference_position == 0 ) continue;
		if ( maximum_number_of_reads_mapped_to_a_single_reference_nucleotide == 0 )
		{
			previous_reference_position = current_reference_position;
			maximum_number_of_reads_mapped_to_a_single_reference_nucleotide = 1;
			number_of_reads_mapped_to_the_current_reference_nucleotide = 1;
		}
		else
		{
			if ( previous_reference_position == current_reference_position )
				number_of_reads_mapped_to_the_current_reference_nucleotide += 1;
			else
			{
				if ( number_of_reads_mapped_to_the_current_reference_nucleotide > maximum_number_of_reads_mapped_to_a_single_reference_nucleotide )
					maximum_number_of_reads_mapped_to_a_single_reference_nucleotide = number_of_reads_mapped_to_the_current_reference_nucleotide;
				previous_reference_position = current_reference_position;
				number_of_reads_mapped_to_the_current_reference_nucleotide = 0;
			}
		}

		if ( strcmp (input_alignment_file_format , "SAM") == 0 )
		{
			line_len = getline ( &line , &len , fhr);
			//printf ("\nLine length %d" , line_len);
		}

		if ( strcmp (input_alignment_file_format , "BAM") == 0 )
			line_len = sam_read1 (fp_in , bamHdr , aln);

		if ( line_len <= 0 ) break;
	} while ( 1 );
	printf ("\nWill start writing to file now");
	fflush (stdout);

	if ( number_of_reads_mapped_to_the_current_reference_nucleotide > maximum_number_of_reads_mapped_to_a_single_reference_nucleotide )
		maximum_number_of_reads_mapped_to_a_single_reference_nucleotide = number_of_reads_mapped_to_the_current_reference_nucleotide;

	strcat(str ,
			"maximum_number_of_reads_mapped_to_a_single_reference_nucleotide: ");
	strcat(str ,
			convertUnsignedIntegerToString (maximum_number_of_reads_mapped_to_a_single_reference_nucleotide));
	strcat(str , "\n");

	strcat(str , "total_number_of_alignments: ");
	strcat(str , convertUnsignedIntegerToString (total_number_of_alignments));
	strcat(str , "\n");

	strcat(str , "max_read_length:");
	strcat(str , convertUnsignedIntegerToString (max_read_length));
	strcat(str , "\n");

	fprintf (fhw , "%s" , str);

	fclose (fhw);
	fclose (fhr);

	if ( strcmp (ended , "BAM") == 0 )
	{
		bam_destroy1 (aln);
		sam_close(fp_in);
	}
}

int main (int argc, char *argv[])
{
	/********************************************************************
	 * Named CLI
	 ********************************************************************/
	struct arguments arguments;

	// Parse our arguments; every option seen by parse_opt will be reflected in arguments.
	// Default values.
	arguments.input_alignment_filename = "";// Empty string - only contains null character
	arguments.summary_information_outputfilename = "";
	arguments.input_alignment_file_format = "";
	arguments.AS_tag_presence = 0;
	arguments.ended = "";

	argp_parse ( &argp , argc , argv , 0 , 0 , &arguments);
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char input_alignment_filename[MAX_FILENAME_LENGTH];
	char input_alignment_file_format[TEN];
	char summary_information_outputfilename[MAX_FILENAME_LENGTH];
	char ended[TEN];
	unsigned short int AS_tag_presence;

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

	strcpy(input_alignment_filename , arguments.input_alignment_filename);
	strcpy(input_alignment_file_format ,
			strupr (arguments.input_alignment_file_format));
	strcpy(summary_information_outputfilename ,
			arguments.summary_information_outputfilename);
	strcpy(ended , arguments.ended);
	AS_tag_presence = arguments.AS_tag_presence;

	/********************************************************************/
	findSummaryInformation (input_alignment_filename ,
			ended ,
			input_alignment_file_format ,
			summary_information_outputfilename ,
			AS_tag_presence);
}
