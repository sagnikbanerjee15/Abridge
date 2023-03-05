#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <argp.h>

// Set up the argument parser
const char *argp_program_version = "abridge printReadNamesAndLineNumbers 1.2.0";
const char *argp_program_bug_address = "sagnikbanerjee15@gmail.com";
static char doc[] = "findAllTags will scan a sam file and output the unique list of alignment tags";
static char args_doc[] = ""; // No standard arguments
							 // (i.e. arguments without "names")

/*
 * Options.  Field 1 in ARGP.
 * Order of fields: {NAME, KEY, ARG, FLAGS, DOC, GROUP}.
 */

static struct argp_option options[] =
	{
		{"inputsamfilename", 'i', "SAMFILENAME", 0, "inputsamfilename of object to download", 0},
		{0, 0, 0, 0, 0, 0} // Last entry should be all zeros in all fields
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
	// char *args[0];   // No standard arguments (without flags)
	char *inputsamfilename; // Argument for --inputsamfilename / -i
};

/*
 * Parser. Field 2 in ARGP.
 * Order of parameters: KEY, ARG, STATE.
 * Parse a single option.
 */

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
	/* Get the input argument from argp_parse, which we
	 know is a pointer to our arguments structure. */
	struct arguments *arguments = state->input;

	// Figure out which option we are parsing, and decide how to store it
	switch (key)
	{
	case 'i':
		arguments->inputsamfilename = arg;
		break;

	case ARGP_KEY_END:
		// Reached the last key.
		// Check if our inputsamfilename and outputfilename REQUIRED "options" have been set to non-default values
		if (strcmp(arguments->inputsamfilename, "") == 0)
		{
			argp_usage(state);
		}
		break;

	default:
		return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

// Our argp parser.
static struct argp argp =
	{options, parse_opt, args_doc, doc, 0, 0, 0};

void printReadNamesAndLineNumbers(char *inputfilename)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	FILE *fhr;
	FILE *fhw;

	size_t len = 0;
	ssize_t line_len;

	char *line = NULL; // for reading each line
	char *read_name = NULL;

	int number_of_fields; // Number of fields in each sam alignment entry

	int i, j, k; // Required in loops

	unsigned long long int line_number;

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	read_name = (char *)malloc(sizeof(char) * 100);

	fhr = fopen(inputfilename, "r");
	if (fhr == NULL)
	{
		printf("Error! File %s not found", inputfilename);
		exit(1);
	}
	/********************************************************************/

	while ((line_len = getline(&line, &len, fhr)) != -1)
		if (line[0] != '@')
			break;
	line_number = 0;
	do
	{
		j = 0;
		for (i = 0; line[i] != '\t'; i++)
			read_name[j++] = line[i];

		printf("%s", read_name);
		printf("\t");
		printf("%lld", line_number);
		printf("\n");
		line_number++;
	} while ((line_len = getline(&line, &len, fhr)) != -1);

	fclose(fhr);
}

int main(int argc, char **argv)
{
	/********************************************************************
	 * Named CLI
	 ********************************************************************/
	struct arguments arguments;

	// Parse our arguments; every option seen by parse_opt will be reflected in arguments.
	// Default values.
	arguments.inputsamfilename = ""; // Empty string - only contains null character

	argp_parse(&argp, argc, argv, 0, 0, &arguments);
	/********************************************************************/
	printReadNamesAndLineNumbers(arguments.inputsamfilename);
	return 0;
}
