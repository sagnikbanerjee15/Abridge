# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <argp.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

// Set up the argument parser
const char *argp_program_version = "abridge findAllTags 1.2.0";
const char *argp_program_bug_address = "sagnikbanerjee15@gmail.com";
static char doc[] = "findAllTags will scan a sam file and output the unique list of alignment tags";
static char args_doc[] = "";  // No standard arguments
							  // (i.e. arguments without "names")

/*
 * Options.  Field 1 in ARGP.
 * Order of fields: {NAME, KEY, ARG, FLAGS, DOC, GROUP}.
 */

static struct argp_option options[] =
{
{ "inputfilename", 'i', "sorted_read_names_with_NH", 0, "Enter the name of the file that has the read names, lines numbers and the NH value", 0 },
{ "outputfilename", 'o', "sorted_read_names_with_NH_and_short_read_names", 0, "Enter the name of the output file. This file will have 4 columns", 0 },
{ "skipshorteningreadname", 's', 0, 0, "Set this option if user requests to skip reducing the size of the read names" },
{ 0, 0, 0, 0, 0, 0 } // Last entry should be all zeros in all fields
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
	// char *args[0];   // No standard arguments (without flags)
	char *inputfilename;	  // Argument for --inputfilename / -i
	char *outputfilename;  // Argument for --outputfilename / -o
	int skipshorteningreadname;
};

/*
 * Parser. Field 2 in ARGP.
 * Order of parameters: KEY, ARG, STATE.
 * Parse a single option.
 */

static error_t parse_opt( int key, char *arg, struct argp_state *state )
{
	/* Get the input argument from argp_parse, which we
	 know is a pointer to our arguments structure. */
	struct arguments *arguments = state->input;

	// Figure out which option we are parsing, and decide how to store it
	switch ( key )
	{
		case 'i':
			arguments->inputfilename = arg;
			break;
		case 'o':
			arguments->outputfilename = arg;
			break;
		case 's':
			arguments->skipshorteningreadname = 1;
			break;

		case ARGP_KEY_END:
			// Reached the last key.
			// Check if our inputsamfilename and outputfilename REQUIRED "options" have been set to non-default values
			if ( strcmp( arguments->inputfilename, "" ) == 0
					|| strcmp( arguments->outputfilename, "" ) == 0 )
			{
				argp_usage( state );
			}
			break;

		default:
			return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

// Our argp parser.
static struct argp argp =
{ options, parse_opt, args_doc, doc, 0, 0, 0 };

void generateNextShortReadId( char *short_read_id )
{
	int number_of_zeros = 0;
	if ( short_read_id[0] == '\0' ) //Empty string
	{
		short_read_id[0] = '0';
		short_read_id[1] = '\0';
	}
	else
	{

		int i;
		i = strlen( short_read_id ) - 1;
		short_read_id[i]++;
		i--;

		for ( ; i >= 0 && short_read_id[i + 1] == 123; i-- )
			short_read_id[i]++;

		for ( i = 0; i < strlen( short_read_id ); i++ )
		{
			if ( short_read_id[i] == 58 )
				short_read_id[i] = 65;
			else if ( short_read_id[i] == 91 )
				short_read_id[i] = 97;
			else if ( short_read_id[i] == 123 )
			{
				short_read_id[i] = 48;
				number_of_zeros++;
			}
		}
		if ( number_of_zeros == strlen( short_read_id ) )
		{
			short_read_id[strlen( short_read_id )] = 48;
			short_read_id[strlen( short_read_id ) + 1] = '\0';
		}
	}
}

void assignShortenedReadsNames(
		char *inputfilename,
		char *outputfilename,
		int skipshorteningreadname )
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;
	FILE *fhw;

	int i;

	size_t len = 0;
	ssize_t line_len;

	char *temp; //Useless
	char *line = NULL; // for reading each line
	char *current_read_id;
	char *previous_read_id;

	char **split_on_tab;
	char *short_read_id;

	unsigned long long int line_number = 0;
	unsigned int number_of_characters_in_short_read_name = 1;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	fhr = fopen( inputfilename, "r" );
	if ( fhr == NULL )
	{
		printf( "\nCannot create file %s", inputfilename );
		exit( 0 );
	}

	fhw = fopen( outputfilename, "w" );
	if ( fhw == NULL )
	{
		printf( "\nCannot create file %s", outputfilename );
		exit( 0 );
	}

	split_on_tab = ( char** ) malloc( sizeof(char*) * 5 );
	for ( i = 0; i < 5; i++ )
		split_on_tab[i] = ( char* ) malloc( sizeof(char) * 1000 );

	current_read_id = ( char* ) malloc( sizeof(char) * 1000 );
	previous_read_id = ( char* ) malloc( sizeof(char) * 1000 );
	short_read_id = ( char* ) malloc( sizeof(char) * 100 );
	short_read_id[0] = '\0';
	current_read_id[0] = '\0';
	previous_read_id[0] = '\0';

	/********************************************************************/

	while ( (line_len = getline( &line, &len, fhr )) != -1 )
	{
		splitByDelimiter( line, '\t', split_on_tab );
		strcpy( current_read_id, split_on_tab[0] );
		if ( skipshorteningreadname == 0 )
		{
			if ( strcmp( previous_read_id, current_read_id ) != 0 ) // Same read name
				generateNextShortReadId( short_read_id );
		}
		else
			strcpy( short_read_id, split_on_tab[0] ); // If shortening of read names is requested to be skipped, then simply copy the long read name

		fprintf( fhw, "%s", split_on_tab[0] );
		fprintf( fhw, "%s", "\t" );
		fprintf( fhw, "%s", split_on_tab[1] );
		fprintf( fhw, "%s", "\t" );
		fprintf( fhw, "%s", split_on_tab[2] );
		fprintf( fhw, "%s", "\t" );
		fprintf( fhw, "%s", short_read_id );
		fprintf( fhw, "%s", "\n" );

		strcpy( previous_read_id, current_read_id );
	}

	fclose( fhr );
	fclose( fhw );
}

int main( int argc, char *argv[] )
{
	/********************************************************************
	 * Named CLI
	 ********************************************************************/
	struct arguments arguments;

	// Parse our arguments; every option seen by parse_opt will be reflected in arguments.
	// Default values.
	arguments.inputfilename = ""; // Empty string - only contains null character
	arguments.outputfilename = "";
	arguments.skipshorteningreadname = 0; // By default read names will be shortened

	argp_parse( &argp, argc, argv, 0, 0, &arguments );
	/********************************************************************/

	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char read_name_to_line_numbers[FILENAME_LENGTH];
	char read_name_to_line_numbers_to_shortened_read_names[FILENAME_LENGTH];

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	strcpy( read_name_to_line_numbers, arguments.inputfilename );
	strcpy(
			read_name_to_line_numbers_to_shortened_read_names,
			arguments.outputfilename );

	/********************************************************************/

	assignShortenedReadsNames(
			read_name_to_line_numbers,
			read_name_to_line_numbers_to_shortened_read_names,
			arguments.skipshorteningreadname );

	return 0;
}
