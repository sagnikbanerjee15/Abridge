#ifndef _ABRIDGE_FUNCTIONS_DEFINITIONS_H_
#define _ABRIDGE_FUNCTIONS_DEFINITIONS_H_

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "data_structure_definitions.h"

/*************************************************************************************************************************
 * START
 *************************************************************************************************************************
 * Functions for allocating memory to Data Structures defined in data_structure_definitions.h
 *************************************************************************************************************************
 *************************************************************************************************************************/

struct Reference_Sequence_Info* allocateMemoryReference_Sequence_Info ()
{
	/*************************************************************************************************************************
	 * Data structure memory allocation function for Reference_Sequence_Info
	 * Initialized memory only for the data structure itself and not for the underlying data
	 * Memory will be dynamically allocated for the data structure variables during file read operation
	 **************************************************************************************************************************/
	struct Reference_Sequence_Info *s;
	s = ( struct Reference_Sequence_Info* ) malloc (sizeof(struct Reference_Sequence_Info));
	s->reference_name = NULL;
	s->nucleotide_sequence = NULL;
	return s;
}

struct Reference* allocateMemoryReference (unsigned int n)
{
	/*************************************************************************************************************************
	 * Data structure memory allocation function for Reference
	 *************************************************************************************************************************/
	struct Reference *s;
	unsigned int i;

	s = ( struct Reference* ) malloc (sizeof(struct Reference));
	s->number_of_reference_sequences = n;
	s->each_sequence = ( struct Reference_Sequence_Info** ) malloc (sizeof(struct Reference_Sequence_Info*) * s->number_of_reference_sequences);
	for ( i = 0 ; i < s->number_of_reference_sequences ; i++ )
	{
		s->each_sequence[i] = allocateMemoryReference_Sequence_Info ();
	}
	return s;
}

struct Cigar_Items* allocateMemoryCigar_Items ()
{
	struct Cigar_Items *instance;
	instance = ( struct Cigar_Items* ) malloc (sizeof(struct Cigar_Items) * 1);
	instance->def = ' ';
	instance->len = 0;
	return instance;
}

struct Samflag_Dictionary_Items* allocateMemorySamflag_Dictionary_Items()
{
	struct Samflag_Dictionary_Items *instance;
	instance = (struct Samflag_Dictionary_Items*) malloc (sizeof(struct Samflag_Dictionary_Items) * 1);
	instance->character = ' ';
	instance->samflag = (char*)malloc(sizeof(char)*TEN);
	return instance;
}

struct Sam_Alignment* allocateMemorySam_Alignment (unsigned int max_read_length)
{
	/********************************************************************
	 * Allocates memory for SAM alignments
	 ********************************************************************/

	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	struct Sam_Alignment *s;
	/********************************************************************/
	s = ( struct Sam_Alignment* ) malloc (sizeof(struct Sam_Alignment));

	max_read_length += 5;
	s->read_name = ( char* ) malloc (sizeof(char) * ONE_HUNDRED);
	s->samflag = ( char* ) malloc (sizeof(char) * TEN);
	s->reference_name = ( char* ) malloc (sizeof(char) * ONE_HUNDRED);
	s->start_position = 0;
	s->mapping_quality_score = ( char* ) malloc (sizeof(char) * TEN);
	s->cigar = ( char* ) malloc (sizeof(char) * max_read_length);
	s->reference_name_next_mate = ( char* ) malloc (sizeof(char) * max_read_length);
	s->start_position_next = 0;
	s->template_length = ( char* ) malloc (sizeof(char) * TEN);
	s->sequence = ( char* ) malloc (sizeof(char) * max_read_length);
	s->quality_scores = ( char* ) malloc (sizeof(char) * max_read_length);
	s->sequence_with_deletions_and_splice_indicators = ( char* ) malloc (sizeof(char) * max_read_length * 2);

	s->level_of_similarity_to_parent_iCIGAR=0;
	s->left_soft_clipped_sequence = ( char* ) malloc (sizeof(char) * max_read_length);
	s->left_soft_clipped_sequence[0] = '\0';
	s->right_soft_clipped_sequence = ( char* ) malloc (sizeof(char) * max_read_length);
	s->right_soft_clipped_sequence[0] = '\0';
	s->left_soft_clipped_qual = ( char* ) malloc (sizeof(char) * max_read_length);
	s->left_soft_clipped_qual[0] = '\0';
	s->right_soft_clipped_qual = ( char* ) malloc (sizeof(char) * max_read_length);
	s->right_soft_clipped_qual[0] = '\0';
	//s->cigar_items = ( struct Cigar_Items* ) malloc (sizeof(struct Cigar_Items) * MAX_CIGAR_ITEMS);
	//s->number_of_cigar_items = 0;

	s->cigar_extended = ( char* ) malloc (sizeof(char) * ( max_read_length * 2 ));
	s->cigar_extended[0] = '\0';
	s->cigar_extended_reference_skips = ( unsigned short int* ) malloc (sizeof(unsigned short int) * ( max_read_length ) * 2);
	s->md_extended = ( char* ) malloc (sizeof(char) * ( max_read_length * 2 ));
	s->md_extended[0] = '\0';
	s->icigar = ( char* ) malloc (sizeof(char) * ( max_read_length * 2 ));
	s->icigar[0] = '\0';
	s->qual_for_mismatches_and_indels = ( char* ) malloc (sizeof(char) * max_read_length);
	s->qual_for_mismatches_and_indels[0] = '\0';
	//s->splices = ( char** ) malloc (sizeof(char*) * 100);
	//for ( i = 0 ; i < 100 ; i++ )
	//	s->splices[i] = ( char* ) malloc (sizeof(char) * 50);
	s->soft_clips_removed_sequence = ( char* ) malloc (sizeof(char) * max_read_length);
	s->soft_clips_removed_sequence[0] = '\0';
	s->soft_clips_removed_quality_scores = ( char* ) malloc (sizeof(char) * max_read_length);
	s->soft_clips_removed_quality_scores[0] = '\0';

	s->NH = ( char* ) malloc (sizeof(char) * 10);
	strcpy(s->NH , "-1");
	s->AS = ( char* ) malloc (sizeof(char) * 5);
	strcpy(s->AS , "X");
	s->MD = ( char* ) malloc (sizeof(char) * ( max_read_length * 2 ));
	strcpy(s->MD , "-1");
	return s;
}

/*************************************************************************************************************************
 *************************************************************************************************************************
 * Functions for allocating memory to Data Structures defined in data_structure_definitions.h
 *************************************************************************************************************************
 * END
 *************************************************************************************************************************/

char* strupr (char *s)
{
	/*************************************************************************************************************************
	 * Convert string to uppercase
	 *************************************************************************************************************************/
	char *pc = s;
	while ( *pc )
	{
		*pc = toupper ( *pc); /* toupper() requires <ctype.h> */
		++pc;
	}
	return s;
}

char* strlwr (char *s)
{
	/*************************************************************************************************************************
	 * Convert string to uppercase
	 *************************************************************************************************************************/
	char *pc = s;
	while ( *pc )
	{
		*pc = tolower ( *pc); /* toupper() requires <ctype.h> */
		++pc;
	}
	return s;
}

void readCompleteReference (
		char *reference_filename,
		struct Reference *reference)
{
	/*************************************************************************************************************************
	 * Opens the reference filename and reports if there are any errors
	 * Operates in 2 passes
	 * 1st pass --> Scans the file to detect the number of references. Reads in only the first character of each line
	 * 2nd pass --> Reads each character in each line and constructs the complete data structure for the reference
	 * Converts lower case nucleotides to upper case
	 *************************************************************************************************************************/

	/*************************************************************************************************************************
	 * Variable declaration and initialization
	 *************************************************************************************************************************/
	FILE *fhr;
	char *buffer = NULL;
	size_t len = 0;
	ssize_t line_len;
	unsigned int i;
	unsigned int j;
	unsigned int number_of_reference_sequences = 0;
	unsigned int number_of_reference_sequences_index = 0;
	/*************************************************************************************************************************/

	fhr = fopen (reference_filename , "rb");
	if ( fhr == NULL )
	{
		printf ("Error! File not found %s" , reference_filename);
		exit (1);
	}

	/*************************************************************************************************************************
	 * 1st pass to determine number of reference sequences
	 *************************************************************************************************************************/
	while ( ( line_len = getline ( &buffer , &len , fhr) ) != -1 )
		if ( buffer[0] == '>' ) number_of_reference_sequences++;

	/*************************************************************************************************************************
	 * Allocate memory for the reference
	 *************************************************************************************************************************/
	reference = allocateMemoryReference (number_of_reference_sequences);

	/*************************************************************************************************************************
	 * 2nd pass for loading all information from the fasta file
	 *************************************************************************************************************************/
	rewind (fhr);
	number_of_reference_sequences_index = -1;

	while ( ( line_len = getline ( &buffer , &len , fhr) ) != -1 )
	{
		// printf("\n%lld", strlen(buffer));
		if ( strlen (buffer) <= 1 ) continue;
		if ( buffer[0] == '>' )
		{
			number_of_reference_sequences_index++;
			reference->each_sequence[number_of_reference_sequences_index]->reference_name = ( char* ) malloc (sizeof(char) * ( line_len + 1 ));
			/*************************************************************************************************************************
			 * Copy the name of the reference only till the first space
			 *************************************************************************************************************************/
			for ( i = 1 , j = 0 ; buffer[i] != 32 ; i++ , j++ )
				reference->each_sequence[number_of_reference_sequences_index]->reference_name[j] = buffer[j];
			reference->each_sequence[number_of_reference_sequences_index]->reference_name[j] = '\0';
		}
		else
		{
			reference->each_sequence[number_of_reference_sequences_index]->nucleotide_sequence = ( char* ) malloc (sizeof(char) * ( line_len + 1 ));
			buffer = strupr (buffer);
			strcpy(reference->each_sequence[number_of_reference_sequences_index]->nucleotide_sequence ,
					buffer);
			reference->each_sequence[number_of_reference_sequences_index]->sequence_length = strlen (buffer);
		}
	}

	fclose (fhr);
}

void convertSignedIntegerToString (char *s, long long int x)
{
	/*************************************************************************************************************************
	 * Converts signed integer to string
	 **************************************************************************************************************************/
	s[0] = '\0';
	sprintf(s , "%lld" , x);
}

void convertUnsignedIntegerToString (char *s, unsigned long long int x)
{
	/*************************************************************************************************************************
	 * Converts unsigned integers to string
	 (*************************************************************************************************************************/
	s[0] = '\0';
	sprintf(s , "%llu" , x);
}

inline unsigned long long int convertStringToUnsignedInteger (char *str)
{
	/*************************************************************************************************************************
	 * Converts string to unsigned integer
	 **************************************************************************************************************************/
	//char convertStringToUnsignedInteger_value[100] = "\0";
	char *convertStringToUnsignedInteger_eptr;

	//convertStringToUnsignedInteger_value[0] = '\0';
	return strtoull (str , &convertStringToUnsignedInteger_eptr , 10);
}

inline long long int convertStringToSignedInteger (char *str)
{
	/*************************************************************************************************************************
	 * Converts string to signed integer
	 **************************************************************************************************************************/
	//char convertStringToSignedInteger_value[100] = "\0";
	char *convertStringToSignedInteger_eptr;

	//convertStringToSignedInteger_value[0] = '\0';
	return strtoll (str , &convertStringToSignedInteger_eptr , 10);
}

unsigned short int splitByDelimiter (
		char *line,
		char delimiter,
		char **new_string)
{
	/*********************************************************************
	 * Splits a string on a requested character and
	 * returns the number of splits
	 *********************************************************************/

	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	int i = 0, j = 0;
	unsigned short ctr = 0;
	/********************************************************************/

	for ( i = 0 ; line[i] != '\0' ; i++ )
	{
		// printf ("\nsplitByDelimiter i=%i" , i);
		//  if space or NULL found, assign NULL into new_string[ctr]
		if ( line[i] == delimiter )
		{
			new_string[ctr][j] = '\0';
			// printf ("\nsplitByDelimiter ctr=%i" , ctr);
			ctr++; // for next word
			j = 0; // for next word, init index to 0
		}
		else if ( line[i] == '\n' )
			continue;
		else
		{
			// printf ("\nsplitByDelimiter j=%i" , j);
			new_string[ctr][j] = line[i];
			j++;
		}
	}
	new_string[ctr][j] = '\0';
	return ctr + 1;
}

int isSequenceSoftClipped (char *cigar)
{
	int i;
	for ( i = 0 ; cigar[i] != '\0' ; i++ )
		if ( cigar[i] == 'S' ) return 1;
	return 0;
}

void splitCigar (char *cigar, // The original CIGAR string
		unsigned short int *num_of_types, //Sets the number of different items in the CIGAR string
		struct Cigar_Items *cigar_items_instance)
{
	/********************************************************************
	 * Variable declarations
	 ********************************************************************/

	int i, j = 0;
	int current_length = 0;
	/********************************************************************/

	for ( i = 0 ; cigar[i] != '\0' ; i++ )
	{
		if ( isdigit (cigar[i]) != 0 ) //cigar[i] is a digit
		{
			current_length = current_length * 10 + cigar[i] - 48;
		}
		else //cigar[i] is a character
		{
			cigar_items_instance[j].def = cigar[i];
			cigar_items_instance[j].len = current_length;
			current_length = 0;
			j += 1;
		}
	}
	*num_of_types = j;
}

void extractSubString (char *str, char *substr, int start_index, int end_index)
{
	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	int length = strlen (str);
	int i;
	int j = 0;

	/********************************************************************/

	if ( length < ( end_index - start_index + 1 ) )
	{
		return;
	}
	for ( i = start_index ; i <= end_index ; i++ )
		substr[j++ ] = str[i];
	substr[j] = '\0';

}

void insertCharacterInString (
		char *str,
		unsigned short int *str_length,
		char ins,
		int loc,
		int number_of_insertions_to_be_made)
{
	/*
	 * Inserts a character. Assumes that a large
	 */
	int i = ( *str_length ) + number_of_insertions_to_be_made;
	for ( ; i > loc ; i-- )
		str[i] = str[i - number_of_insertions_to_be_made];
	( *str_length ) += number_of_insertions_to_be_made;
	while ( number_of_insertions_to_be_made-- )
		str[i++ ] = ins;
}

char returnEncodedInsertionCharacter (char nucleotide)
{
	switch ( nucleotide )
	{
		case 'A':
			return '!';
		case 'T':
			return '"';
		case 'G':
			return '#';
		case 'C':
			return '$';
		case 'N':
			return '%';
	}
	return ' ';
}

char returnEncodedMismatchCharacter (char nucleotide)
{
	switch ( nucleotide )
	{
		case 'A':
			return '&';
		case 'T':
			return '\'';
		case 'G':
			return '(';
		case 'C':
			return ')';
		case 'N':
			return '*';
	}
	return ' ';
}

void generateiCIGARString (
		struct Sam_Alignment *sam_alignment_instance,
		unsigned short int AS_tag_presence,
		unsigned short int flag_ignore_alignment_scores,
		unsigned short int flag_ignore_soft_clippings,
		unsigned short int flag_ignore_mismatches,
		unsigned short int flag_ignore_all_quality_scores,
		unsigned short int flag_ignore_unmapped_sequences,
		unsigned short int flag_ignore_quality_scores_for_matched_bases,
		char *ended,
		struct Cigar_Items *cigar_items_instance,
		struct Samflag_Dictionary_Items **samflag_dictionary,
		unsigned int samflag_dictionary_size)
{
	/*************************************************************************************************************************
	 * Generate the integrated CIGAR string from CIGAR and MD
	 *************************************************************************************************************************/

	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	int left_soft_clip_point = 0;
	int right_soft_clip_point = 0;
	int flag;
	int print_outputs = 0;
	int perfect_alignment_indicator = 0;
	//int spliced_alignment_indicator = 0;
	char str[1000];
	char temp_str[5];
	char M_replacement_character;
	char dummy;

	unsigned short int cigar_items_index;
	unsigned short int total_number_of_cigar_items;
	unsigned short int expanded_cigar_string_index;
	unsigned short int expanded_cigar_string_length;
	unsigned short int insertion_nucleotides_start_index;
	unsigned short int sequence_with_deletions_and_splice_indicators_length;
	unsigned short int index_for_N_and_D_insertion_in_sequence;

	unsigned short int cigar_expansion_iterator = 0;
	unsigned short int MD_index = 0;

	unsigned short int expanded_md_string_index = 0;
	unsigned short int expanded_md_string_length = 0;

	unsigned short int icigar_index;
	unsigned short int icigar_length;

	unsigned short int num;
	unsigned short int i, j;

	unsigned short int runlength;
	unsigned int samflag_dictionary_index;
	/********************************************************************/
	sequence_with_deletions_and_splice_indicators_length = strlen (sam_alignment_instance->sequence_with_deletions_and_splice_indicators);
	expanded_cigar_string_index = 0;
	sam_alignment_instance->cigar_extended[0] = '\0';

	splitCigar (sam_alignment_instance->cigar ,
			&total_number_of_cigar_items ,
			cigar_items_instance);

	insertion_nucleotides_start_index = 0;
	index_for_N_and_D_insertion_in_sequence = 0;
	for ( cigar_items_index = 0 ;
			cigar_items_index < total_number_of_cigar_items ;
			cigar_items_index++ )
	{
		switch ( cigar_items_instance[cigar_items_index].def )
		{
			case 'S':
				for ( cigar_expansion_iterator = 0 ;
						cigar_expansion_iterator < cigar_items_instance[cigar_items_index].len ;
						cigar_expansion_iterator++ )
				{
					sam_alignment_instance->cigar_extended[expanded_cigar_string_index] = 'S';
					sam_alignment_instance->cigar_extended_reference_skips[expanded_cigar_string_index] = 1;
					expanded_cigar_string_index++;
				}
				insertion_nucleotides_start_index += cigar_items_instance[cigar_items_index].len;
				index_for_N_and_D_insertion_in_sequence += cigar_items_instance[cigar_items_index].len;
				break;
			case 'M':
				for ( cigar_expansion_iterator = 0 ;
						cigar_expansion_iterator < cigar_items_instance[cigar_items_index].len ;
						cigar_expansion_iterator++ )
				{
					sam_alignment_instance->cigar_extended[expanded_cigar_string_index] = 'M';
					sam_alignment_instance->cigar_extended_reference_skips[expanded_cigar_string_index] = 1;
					expanded_cigar_string_index++;
				}
				insertion_nucleotides_start_index += cigar_items_instance[cigar_items_index].len;
				index_for_N_and_D_insertion_in_sequence += cigar_items_instance[cigar_items_index].len;
				break;
			case '=':
				for ( cigar_expansion_iterator = 0 ;
						cigar_expansion_iterator < cigar_items_instance[cigar_items_index].len ;
						cigar_expansion_iterator++ )
				{
					sam_alignment_instance->cigar_extended[expanded_cigar_string_index] = 'M';
					sam_alignment_instance->cigar_extended_reference_skips[expanded_cigar_string_index] = 1;
					expanded_cigar_string_index++;
				}
				insertion_nucleotides_start_index += cigar_items_instance[cigar_items_index].len;
				index_for_N_and_D_insertion_in_sequence += cigar_items_instance[cigar_items_index].len;
				break;
			case 'X':
				for ( cigar_expansion_iterator = 0 ;
						cigar_expansion_iterator < cigar_items_instance[cigar_items_index].len ;
						cigar_expansion_iterator++ )
				{
					sam_alignment_instance->cigar_extended[expanded_cigar_string_index++ ] = 'M';
					sam_alignment_instance->cigar_extended_reference_skips[expanded_cigar_string_index] = 1;
					expanded_cigar_string_index++;
				}
				insertion_nucleotides_start_index += cigar_items_instance[cigar_items_index].len;
				index_for_N_and_D_insertion_in_sequence += cigar_items_instance[cigar_items_index].len;
				break;
			case 'I':
				for ( cigar_expansion_iterator = 0 ;
						cigar_expansion_iterator < cigar_items_instance[cigar_items_index].len ;
						cigar_expansion_iterator++ )
				{
					sam_alignment_instance->cigar_extended[expanded_cigar_string_index] = returnEncodedInsertionCharacter (sam_alignment_instance->sequence[insertion_nucleotides_start_index + cigar_expansion_iterator]);
					sam_alignment_instance->cigar_extended_reference_skips[expanded_cigar_string_index] = 1;
					expanded_cigar_string_index++;
				}
				insertion_nucleotides_start_index += cigar_items_instance[cigar_items_index].len;
				index_for_N_and_D_insertion_in_sequence += cigar_items_instance[cigar_items_index].len;
				break;
			case 'D':
			{
				sam_alignment_instance->cigar_extended[expanded_cigar_string_index] = 'D';
				sam_alignment_instance->cigar_extended_reference_skips[expanded_cigar_string_index] = cigar_items_instance[cigar_items_index].len;
				expanded_cigar_string_index++;

				insertCharacterInString (sam_alignment_instance->sequence_with_deletions_and_splice_indicators ,
						&sequence_with_deletions_and_splice_indicators_length ,
						'D' ,
						index_for_N_and_D_insertion_in_sequence ,
						1);
			}
				index_for_N_and_D_insertion_in_sequence += 1;
				break;
			case 'N':
			{
				sam_alignment_instance->cigar_extended[expanded_cigar_string_index] = 'N';
				sam_alignment_instance->cigar_extended_reference_skips[expanded_cigar_string_index] = cigar_items_instance[cigar_items_index].len;
				expanded_cigar_string_index++;

				insertCharacterInString (sam_alignment_instance->sequence_with_deletions_and_splice_indicators ,
						&sequence_with_deletions_and_splice_indicators_length ,
						'X' ,
						index_for_N_and_D_insertion_in_sequence ,
						1);
			}
				index_for_N_and_D_insertion_in_sequence += 1;
				break;
		}
	}
	expanded_cigar_string_length = expanded_cigar_string_index;
	sam_alignment_instance->cigar_extended[expanded_cigar_string_length] = '\0';

	/************************************************************************************************************************
	 * Expand the MD string
	 ************************************************************************************************************************/
	MD_index = 0;
	num = 0;

	// Check if left soft clip exists
	if ( cigar_items_instance[0].def == 'S' )
	{
		for ( cigar_expansion_iterator = 0 ;
				cigar_expansion_iterator < cigar_items_instance[0].len ;
				cigar_expansion_iterator++ )
		{
			sam_alignment_instance->md_extended[expanded_md_string_index++ ] = 'S';
		}
	}
	while ( sam_alignment_instance->MD[MD_index] != '\0' )
	{
		if ( isdigit (sam_alignment_instance->MD[MD_index]) != 0 ) // md[i] is a digit
			num = num * 10 + sam_alignment_instance->MD[MD_index] - 48;
		else if ( sam_alignment_instance->MD[MD_index] == '^' )
		{
			for ( j = 0 ; j < num ; j++ )
				sam_alignment_instance->md_extended[expanded_md_string_index++ ] = '=';
			num = 0;
			MD_index += 1;
			while ( isdigit (sam_alignment_instance->MD[MD_index]) == 0 ) // Iterate till you pick up all the deletions and ignore them
				MD_index += 1;
			MD_index -= 1;
		}
		else
		{
			for ( j = 0 ; j < num ; j++ )
				sam_alignment_instance->md_extended[expanded_md_string_index++ ] = '=';
			sam_alignment_instance->md_extended[expanded_md_string_index++ ] = returnEncodedMismatchCharacter (sam_alignment_instance->MD[MD_index]);
			num = 0;
		}
		MD_index += 1;
	}
	for ( j = 0 ; j < num ; j++ )
		sam_alignment_instance->md_extended[expanded_md_string_index++ ] = '=';

	// Check if right soft clip exists
	if ( cigar_items_instance[total_number_of_cigar_items - 1].def == 'S' )
	{
		for ( cigar_expansion_iterator = 0 ;
				cigar_expansion_iterator < cigar_items_instance[total_number_of_cigar_items - 1].len ;
				cigar_expansion_iterator++ )
		{
			sam_alignment_instance->md_extended[expanded_md_string_index++ ] = 'S';
		}
	}
	sam_alignment_instance->md_extended[expanded_md_string_index] = '\0';
	expanded_md_string_length = expanded_md_string_index;

	/*************************************************************************************************************************
	 * Adjust the expanded MD string to accommodate insert characters, N and D
	 *************************************************************************************************************************/
	j = 0;
	for ( i = 0 ; i < expanded_cigar_string_length ; i++ )
	{
		if ( sam_alignment_instance->cigar_extended[i] >= 33 && sam_alignment_instance->cigar_extended[i] <= 37 )
		{
			insertCharacterInString (sam_alignment_instance->md_extended ,
					&expanded_md_string_length ,
					sam_alignment_instance->cigar_extended[i] ,
					i ,
					1);
			j = -1;
		}
		else if ( sam_alignment_instance->cigar_extended[i] == 'N' || sam_alignment_instance->cigar_extended[i] == 'D' )
		{
			insertCharacterInString (sam_alignment_instance->md_extended ,
					&expanded_md_string_length ,
					sam_alignment_instance->cigar_extended[i] ,
					i ,
					1);
		}
	}
	/************************************************************************************************************************/

	/*************************************************************************************************************************
	 * Make final adjustments to insert and mismatch characters to represent the nucleotide in the read and not in the reference
	 *************************************************************************************************************************/
	for ( i = 0 ; i < expanded_cigar_string_length ; i++ )
	{
		switch ( sam_alignment_instance->md_extended[i] )
		{
			//Insert characters:
			case '!':
			case '"':
			case '#':
			case '$':
			case '%':
				sam_alignment_instance->md_extended[i] = returnEncodedInsertionCharacter (sam_alignment_instance->sequence_with_deletions_and_splice_indicators[i]);
				break;

				//Mismatch characters
			case '&':
			case '\'':
			case '(':
			case ')':
			case '*':
				sam_alignment_instance->md_extended[i] = returnEncodedMismatchCharacter (sam_alignment_instance->sequence_with_deletions_and_splice_indicators[i]);
				break;

		}
	}
	/************************************************************************************************************************/

	/*************************************************************************************************************************
	 * Construct the iCIGAR from all the information
	 *************************************************************************************************************************/
	icigar_index = 0;
	runlength = 0;
	sam_alignment_instance->icigar[0] = '\0';
	for ( i = 0 ; i < expanded_cigar_string_length ; i++ )
	{

		if ( i == 0 && sam_alignment_instance->cigar_extended[i] == 'S' ) //Check for left soft clip
		{
			if ( flag_ignore_soft_clippings == 0 )
			{
				strcat(sam_alignment_instance->icigar ,
						strlwr (sam_alignment_instance->left_soft_clipped_sequence));
			}
			while ( sam_alignment_instance->cigar_extended[i] == 'S' )
				i++;
			i--;
		}
		else if ( sam_alignment_instance->cigar_extended[i] == 'S' ) //Check for right soft clip
		{
			if ( runlength != 0 )
			{
				convertUnsignedIntegerToString (str ,
						( unsigned long long ) runlength);
				strcat(sam_alignment_instance->icigar , str);

				runlength = 0;

				str[0] = 'M';
				str[1] = '\0';
				strcat(sam_alignment_instance->icigar , str);
			}
			if ( flag_ignore_soft_clippings == 0 )
			{
				strcat(sam_alignment_instance->icigar ,
						strlwr (sam_alignment_instance->right_soft_clipped_sequence));
			}
			break;
		}
		else if ( sam_alignment_instance->cigar_extended[i] == 'M' ) // Match or mismatch
		{
			if ( sam_alignment_instance->md_extended[i] == '=' )
				runlength++; // a match
			else // mismatch found
			{
				if ( runlength != 0 )
				{
					convertUnsignedIntegerToString (str ,
							( unsigned long long ) runlength);
					strcat(sam_alignment_instance->icigar , str);

					runlength = 0;

					str[0] = 'M';
					str[1] = '\0';
					strcat(sam_alignment_instance->icigar , str);
				}

				str[0] = sam_alignment_instance->md_extended[i];
				str[1] = '\0';
				strcat(sam_alignment_instance->icigar , str);

			}
		}
		else if ( sam_alignment_instance->cigar_extended[i] == 'D' ) // Deletion from reference
		{
			if ( runlength != 0 )
			{
				convertUnsignedIntegerToString (str ,
						( unsigned long long ) runlength);
				strcat(sam_alignment_instance->icigar , str);

				runlength = 0;

				str[0] = 'M';
				str[1] = '\0';
				strcat(sam_alignment_instance->icigar , str);
			}

			convertUnsignedIntegerToString (str ,
					( unsigned long long ) sam_alignment_instance->cigar_extended_reference_skips[i]);
			strcat(sam_alignment_instance->icigar , str);

			str[0] = 'D';
			str[1] = '\0';
			strcat(sam_alignment_instance->icigar , str);
		}
		else if ( sam_alignment_instance->cigar_extended[i] == 'N' ) // Splices
		{
			if ( runlength != 0 )
			{
				convertUnsignedIntegerToString (str ,
						( unsigned long long ) runlength);
				strcat(sam_alignment_instance->icigar , str);

				runlength = 0;

				str[0] = 'M';
				str[1] = '\0';
				strcat(sam_alignment_instance->icigar , str);
			}
			convertUnsignedIntegerToString (str ,
					( unsigned long long ) sam_alignment_instance->cigar_extended_reference_skips[i]);
			strcat(sam_alignment_instance->icigar , str);

			str[0] = 'N';
			str[1] = '\0';
			strcat(sam_alignment_instance->icigar , str);
		}
		else // Insertions
		{
			if ( runlength != 0 )
			{
				convertUnsignedIntegerToString (str ,
						( unsigned long long ) runlength);
				strcat(sam_alignment_instance->icigar , str);

				runlength = 0;

				str[0] = 'M';
				str[1] = '\0';
				strcat(sam_alignment_instance->icigar , str);
			}

			str[0] = sam_alignment_instance->cigar_extended[i];
			str[1] = '\0';
			strcat(sam_alignment_instance->icigar , str);
		}
	}

	if ( runlength != 0 )
	{
		convertUnsignedIntegerToString (str , ( unsigned long long ) runlength);
		strcat(sam_alignment_instance->icigar , str);

		runlength = 0;

		str[0] = 'M';
		str[1] = '\0';
		strcat(sam_alignment_instance->icigar , str);
	}

	// Replace the M in the iCIGAR with samformatflag character

	for(samflag_dictionary_index=0; samflag_dictionary_index<samflag_dictionary_size;samflag_dictionary_index++)
	{
		//printf("\nsamflag_in_dictionary=%s samflag_in_alignment=%s replacement_character = %c",samflag_dictionary[samflag_dictionary_index]->samflag, sam_alignment_instance->samflag, samflag_dictionary[samflag_dictionary_index]->character);
		if(strcmp(samflag_dictionary[samflag_dictionary_index]->samflag, sam_alignment_instance->samflag) == 0)
			break;
	}
	if(samflag_dictionary_index==samflag_dictionary_size)
	{
		printf("\nBig Trouble");
	}
	sam_alignment_instance->replacement_character = samflag_dictionary[samflag_dictionary_index]->character;

	//Add the NH tag, NH is stored as a string
	strcpy(sam_alignment_instance->icigar,sam_alignment_instance->NH);

	// Add the MAPQ and the AS scores
	if(flag_ignore_alignment_scores == 0)
	{
		strcpy(sam_alignment_instance->icigar,"~");
		strcpy(sam_alignment_instance->icigar,sam_alignment_instance->mapping_quality_score);
		strcpy(sam_alignment_instance->icigar,"~");
		strcpy(sam_alignment_instance->icigar,sam_alignment_instance->AS);
	}
	/************************************************************************************************************************/

	/*printf ("\nCIGAR=%s\tMD=%s\n%s\n%s\n%s\n%s" ,
			sam_alignment_instance->cigar ,
			sam_alignment_instance->MD ,
			sam_alignment_instance->cigar_extended ,
			sam_alignment_instance->md_extended ,
			sam_alignment_instance->sequence_with_deletions_and_splice_indicators ,
			sam_alignment_instance->icigar);
	fflush (stdout);
	*/

}

void replaceSingleCharacterInString(char *str, char old_char, char new_char)
{
	if(old==new)
		return;
	for(unsigned int i=0;str[i]!='\0';i++)
	{
		if(str[i]==old_char)
			str[i]=new_char;
	}
}

void processSoftClips (
		struct Sam_Alignment *dest,
		unsigned short int AS_tag_presence,
		unsigned short int flag_ignore_alignment_scores,
		unsigned short int flag_ignore_soft_clippings,
		unsigned short int flag_ignore_mismatches,
		unsigned short int flag_ignore_all_quality_scores,
		unsigned short int flag_ignore_unmapped_sequences,
		unsigned short int flag_ignore_quality_scores_for_matched_bases,
		unsigned int max_read_length,
		struct Cigar_Items *cigar_items_instance)
{

	unsigned short int total_number_of_cigar_items;
	unsigned short int left_soft_clip_point;
	unsigned short int right_soft_clip_point;

	unsigned int i;
	/*************************************************************************************************************************
	 * Processes the Sam_Alignemnt data structure and constructs the soft clips
	 *************************************************************************************************************************/
	if ( isSequenceSoftClipped (dest->cigar) == 0 || flag_ignore_soft_clippings == 1 ) // Process everything when there are no soft clips
	{
		strcpy(dest->left_soft_clipped_sequence , "");
		strcpy(dest->right_soft_clipped_sequence , "");
		strcpy(dest->left_soft_clipped_qual , "");
		strcpy(dest->right_soft_clipped_qual , "");
		dest->left_soft_clipped_sequence_length = 0;
		dest->right_soft_clipped_sequence_length = 0;
		strcpy(dest->soft_clips_removed_sequence , dest->sequence);
		strcpy(dest->soft_clips_removed_quality_scores , dest->quality_scores);
		dest->soft_clips_removed_sequence_len = strlen (dest->soft_clips_removed_sequence);

	}
	else
	{
		splitCigar (dest->cigar ,
				&total_number_of_cigar_items ,
				cigar_items_instance);

		left_soft_clip_point = 0;
		right_soft_clip_point = 0;
		if ( cigar_items_instance[0].def == 'S' )
		{
			left_soft_clip_point = cigar_items_instance[0].len;
			extractSubString (dest->sequence ,
					dest->left_soft_clipped_sequence ,
					0 ,
					left_soft_clip_point - 1);
			extractSubString (dest->quality_scores ,
					dest->left_soft_clipped_qual ,
					0 ,
					left_soft_clip_point - 1);
			dest->left_soft_clipped_sequence_length = dest->read_sequence_len - cigar_items_instance[0].len;

		}
		if ( cigar_items_instance[total_number_of_cigar_items - 1].def == 'S' )
		{
			right_soft_clip_point = cigar_items_instance[total_number_of_cigar_items - 1].len;
			extractSubString (dest->sequence ,
					dest->right_soft_clipped_sequence ,
					dest->read_sequence_len - right_soft_clip_point ,
					dest->read_sequence_len - 1);
			extractSubString (dest->quality_scores ,
					dest->right_soft_clipped_qual ,
					dest->read_sequence_len - right_soft_clip_point ,
					dest->read_sequence_len - 1);
			dest->right_soft_clipped_sequence_length = dest->read_sequence_len - cigar_items_instance[total_number_of_cigar_items - 1].len;
		}

		/*
		 * Remove the soft-clips and construct new sequence
		 */

		int j = 0;
		for ( i = left_soft_clip_point ;
				i < dest->read_sequence_len - right_soft_clip_point ; i++ )
		{
			dest->soft_clips_removed_sequence[j] = dest->sequence[i];
			dest->soft_clips_removed_quality_scores[j] = dest->quality_scores[i];
			j++;
		}
		dest->soft_clips_removed_sequence[j] = '\0';
		dest->soft_clips_removed_quality_scores[j] = '\0';
		dest->soft_clips_removed_sequence_len = strlen (dest->soft_clips_removed_sequence);
	}
}

void populateSamAlignmentInstance (
		struct Sam_Alignment *dest,
		char **src,
		char **split_on_colon,
		int number_of_fields,
		unsigned short int AS_tag_presence,
		unsigned short int flag_ignore_alignment_scores,
		unsigned short int flag_ignore_soft_clippings,
		unsigned short int flag_ignore_mismatches,
		unsigned short int flag_ignore_all_quality_scores,
		unsigned short int flag_ignore_unmapped_sequences,
		unsigned short int flag_ignore_quality_scores_for_matched_bases)
{
	/*************************************************************************************************************************
	 * Processed the line read from the SAM alignment file and prepares the Sam_Alignment instance
	 * Skips unimportant fields if requested
	 *************************************************************************************************************************/

	/********************************************************************
	 * Variable declarations
	 ********************************************************************/
	char *temp;
	int i;

	/********************************************************************/

	strcpy(dest->read_name , src[0]);
	strcpy(dest->samflag, src[1]);
	strcpy(dest->reference_name , src[2]);
	dest->start_position = convertStringToUnsignedInteger (src[3]);
	strcpy(dest->mapping_quality_score , src[4]);
	strcpy(dest->cigar , src[5]);
	strcpy(dest->reference_name_next_mate , src[6]);
	dest->start_position_next = convertStringToUnsignedInteger (src[7]);
	strcpy(dest->template_length , src[8]);
	strcpy(dest->sequence , strupr (src[9]));
	strcpy(dest->sequence_with_deletions_and_splice_indicators ,
			dest->sequence);
	strcpy(dest->quality_scores , src[10]);
	//for ( i = 0 ; dest->quality_scores[i] != '\0' ; i++ )
	//	dest->quality_scores[i] += QUAL_SCORE_ADJUSTMENT;

	for ( i = 11 ; i < number_of_fields ; i++ )
	{
		//printf("%s\n",src[i]);
		//fflush ( stdout );
		splitByDelimiter (src[i] , ':' , split_on_colon);
		//sam_tag_index = locateSamTags(split_on_colon[0]);
		//dest->tags[i - 11].name = sam_tags[sam_tag_index];
		if ( i == number_of_fields - 1 )
			split_on_colon[2][strlen (split_on_colon[2])] = '\0';

		if ( strcmp (split_on_colon[0] , "NH") == 0 )
			strcpy(dest->NH , split_on_colon[2]);
		else if ( strcmp (split_on_colon[0] , "MD") == 0 )
			strcpy(dest->MD , split_on_colon[2]);
		else if ( strcmp (split_on_colon[0] , "AS") == 0 )
		{
			strcpy(dest->AS , split_on_colon[2]);
		}

		//printf("\n Tags %s Parts of the tag %s %s %s ", src[i], split_on_colon[0], split_on_colon[1], split_on_colon[2]);
	}

	dest->read_sequence_len = strlen (dest->sequence);

}

void populateBamAlignmentInstance (struct Sam_Alignment *dest, // Final Sam alignment DS
		samFile *fp_in, // File pointer if BAM file provided
		bam_hdr_t *bamHdr,         // read header
		bam1_t *aln,
		char **split_on_colon,
		int number_of_fields,
		unsigned short int AS_tag_presence,
		unsigned short int flag_ignore_alignment_scores,
		unsigned short int flag_ignore_soft_clippings,
		unsigned short int flag_ignore_mismatches,
		unsigned short int flag_ignore_all_quality_scores,
		unsigned short int flag_ignore_unmapped_sequences,
		unsigned short int flag_ignore_quality_scores_for_matched_bases)
{

	int i;
	strcpy(dest->read_name , bam_get_qname (aln));
	dest->samflag = aln->core.flag;
	strcpy(dest->reference_name , bamHdr->target_name[aln->core.tid]); // contig name (chromosome)
	dest->start_position = aln->core.pos + 1; // left most position of alignment in zero based coordianate (+1)
	dest->mapping_quality_score = aln->core.qual;

	uint32_t len = aln->core.l_qseq;                // length of the read.

	uint8_t *q = bam_get_seq (aln); // quality string

	char *qseq = ( char* ) malloc (len);

	for ( i = 0 ; i < len ; i++ )
	{
		qseq[i] = seq_nt16_str[bam_seqi (q , i)]; // gets nucleotide id and converts them into IUPAC id.
	}
	qseq[i] = '\0';
}

unsigned short int prepareSingleRecordFromAlignmentFile (
		char *line,
		samFile *fp_in,           // File pointer if BAM file provided
		bam_hdr_t *bamHdr,           // read header
		bam1_t *aln,
		FILE *fhr,
		struct Sam_Alignment *sam_alignment_instance,
		char *ended,
		char *alignment_format,
		unsigned short int AS_tag_presence,
		unsigned short int flag_ignore_alignment_scores,
		unsigned short int flag_ignore_soft_clippings,
		unsigned short int flag_ignore_mismatches,
		unsigned short int flag_ignore_all_quality_scores,
		unsigned short int flag_ignore_unmapped_sequences,
		unsigned short int flag_ignore_quality_scores_for_matched_bases,
		unsigned int max_read_length,
		char **split_on_tab,
		char **split_on_colon,
		struct Cigar_Items *cigar_items_instance,
		struct Samflag_Dictionary_Items **samflag_dictionary,
		unsigned int samflag_dictionary_size)
{
	/*************************************************************************************************************************
	 * Reads in each alignment and stores the values in them Sam_Alignment object
	 * Decided which fields to retain based on whether the alignment was SE or PE
	 **************************************************************************************************************************/

	size_t len = 0;
	ssize_t line_len;

	unsigned short int number_of_fields;

	if ( strcmp (alignment_format , "SAM") == 0 )
	{
		//line_len = getline ( &line , &len , fhr);
		//if ( line_len == -1 ) return -1;
		number_of_fields = splitByDelimiter (line , '\t' , split_on_tab);
		populateSamAlignmentInstance (sam_alignment_instance ,
				split_on_tab ,
				split_on_colon ,
				number_of_fields ,
				AS_tag_presence ,
				flag_ignore_alignment_scores ,
				flag_ignore_soft_clippings ,
				flag_ignore_mismatches ,
				flag_ignore_all_quality_scores ,
				flag_ignore_unmapped_sequences ,
				flag_ignore_quality_scores_for_matched_bases);
	}
	else if ( strcmp (alignment_format , "BAM") == 0 )
	{
		populateBamAlignmentInstance (sam_alignment_instance , fp_in ,// File pointer if BAM file provided
				bamHdr ,		// read header
				aln ,
				split_on_colon ,
				number_of_fields ,
				AS_tag_presence ,
				flag_ignore_alignment_scores ,
				flag_ignore_soft_clippings ,
				flag_ignore_mismatches ,
				flag_ignore_all_quality_scores ,
				flag_ignore_unmapped_sequences ,
				flag_ignore_quality_scores_for_matched_bases);
	}
	/*
	 * Process soft clipped fields - find the soft clipped portions of the reads and prepare the corresponding fields
	 */
	processSoftClips (sam_alignment_instance ,
			AS_tag_presence ,
			flag_ignore_alignment_scores ,
			flag_ignore_soft_clippings ,
			flag_ignore_mismatches ,
			flag_ignore_all_quality_scores ,
			flag_ignore_unmapped_sequences ,
			flag_ignore_quality_scores_for_matched_bases ,
			max_read_length ,
			cigar_items_instance);

	/*
	 * Generate iCIGAR string
	 */
	generateiCIGARString (sam_alignment_instance ,
			AS_tag_presence ,
			flag_ignore_alignment_scores ,
			flag_ignore_soft_clippings ,
			flag_ignore_mismatches ,
			flag_ignore_all_quality_scores ,
			flag_ignore_unmapped_sequences ,
			flag_ignore_quality_scores_for_matched_bases ,
			ended ,
			cigar_items_instance,
			samflag_dictionary,
			samflag_dictionary_size);
	return 0;
}

#endif /* ABRIDGE_FUNCTIONS_DEFINITIONS_H_ */
