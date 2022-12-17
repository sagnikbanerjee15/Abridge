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

struct Sam_Alignment* allocateMemorySam_Alignment ()
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

	s->read_name = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->reference_name = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->start_position = 0;
	s->mapping_quality_score = ( char* ) malloc (sizeof(char) * TEN);
	s->cigar = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->reference_name_next_mate = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->start_position_next = 0;
	s->template_length = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->sequence = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->quality_scores = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);

	s->left_soft_clipped_sequence = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->left_soft_clipped_sequence[0] = '\0';
	s->right_soft_clipped_sequence = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->right_soft_clipped_sequence[0] = '\0';
	s->left_soft_clipped_qual = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->left_soft_clipped_qual[0] = '\0';
	s->right_soft_clipped_qual = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->right_soft_clipped_qual[0] = '\0';
	//s->cigar_items = ( struct Cigar_Items* ) malloc (sizeof(struct Cigar_Items) * MAX_CIGAR_ITEMS);
	//s->number_of_cigar_items = 0;

	s->cigar_extended = ( char* ) malloc (sizeof(char) * ( MAX_SEQ_LEN * 2 ));
	s->cigar_extended[0] = '\0';
	s->md_extended = ( char* ) malloc (sizeof(char) * ( MAX_SEQ_LEN * 2 ));
	s->md_extended[0] = '\0';
	s->icigar = ( char* ) malloc (sizeof(char) * ( MAX_SEQ_LEN * 2 ));
	s->icigar[0] = '\0';
	//s->splices = ( char** ) malloc (sizeof(char*) * 100);
	//for ( i = 0 ; i < 100 ; i++ )
	//	s->splices[i] = ( char* ) malloc (sizeof(char) * 50);
	s->soft_clips_removed_sequence = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->soft_clips_removed_sequence[0] = '\0';
	s->soft_clips_removed_quality_scores = ( char* ) malloc (sizeof(char) * MAX_SEQ_LEN);
	s->soft_clips_removed_quality_scores[0] = '\0';

	s->NH = ( char* ) malloc (sizeof(char) * 10);
	strcpy(s->NH , "-1");
	s->AS = ( char* ) malloc (sizeof(char) * 5);
	strcpy(s->AS , "X");
	s->MD = ( char* ) malloc (sizeof(char) * ( MAX_SEQ_LEN * 2 ));
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

inline char* convertSignedIntegerToString (long long int x)
{
	/*************************************************************************************************************************
	 * Converts signed integer to string
	 **************************************************************************************************************************/
	char str[256];
	sprintf(str , "%lld" , x);

	return str;
}

inline char* convertUnsignedIntegerToString (unsigned long long int x)
{
	/*************************************************************************************************************************
	 * Converts unsigned integers to string
	 (*************************************************************************************************************************/
	char str[256];
	sprintf(str , "%llu" , x);

	return str;
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

void populateSamAlignmentInstance (
		struct Sam_Alignment *dest,
		char **src,
		char **split_tags,
		int number_of_fields,
		unsigned short int AS_tag_presence)
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

	struct Cigar_Items cigar_items_instance[MAX_SEQ_LEN];
	unsigned short int total_number_of_cigar_items;
	unsigned short int left_soft_clip_point;
	unsigned short int right_soft_clip_point;
	/********************************************************************/

	strcpy(dest->read_name , src[0]);
	dest->samflag = convertStringToUnsignedInteger (src[1]);
	strcpy(dest->reference_name , src[2]);
	dest->start_position = convertStringToUnsignedInteger (src[3]);
	strcpy(dest->mapping_quality_score , src[4]);
	strcpy(dest->cigar , src[5]);
	strcpy(dest->reference_name_next_mate , src[6]);
	dest->start_position_next = convertStringToUnsignedInteger (src[7]);
	strcpy(dest->template_length , src[8]);
	strcpy(dest->sequence , strupr (src[9]));
	strcpy(dest->quality_scores , src[10]);
	for ( i = 0 ; dest->quality_scores[i] != '\0' ; i++ )
		dest->quality_scores[i] += QUAL_SCORE_ADJUSTMENT;

	for ( i = 11 ; i < number_of_fields ; i++ )
	{
		//printf("%s\n",src[i]);
		//fflush ( stdout );
		splitByDelimiter (src[i] , ':' , split_tags);
		//sam_tag_index = locateSamTags(split_tags[0]);
		//dest->tags[i - 11].name = sam_tags[sam_tag_index];
		if ( i == number_of_fields - 1 )
			split_tags[2][strlen (split_tags[2])] = '\0';

		if ( strcmp (split_tags[0] , "NH") == 0 )
			strcpy(dest->NH , split_tags[2]);
		else if ( strcmp (split_tags[0] , "MD") == 0 )
			strcpy(dest->MD , split_tags[2]);
		else if ( strcmp (split_tags[0] , "AS") == 0 )
		{
			strcpy(dest->AS , split_tags[2]);
		}

		//printf("\n Tags %s Parts of the tag %s %s %s ", src[i], split_tags[0], split_tags[1], split_tags[2]);
	}

	dest->read_sequence_len = strlen (dest->sequence);

	/*
	 * Process soft clipped fields
	 */
	if ( isSequenceSoftClipped (dest->cigar) == 0 ) // Process everything when there are no soft clips
	{
		strcpy(dest->left_soft_clipped_sequence , "");
		strcpy(dest->right_soft_clipped_sequence , "");
		strcpy(dest->left_soft_clipped_qual , "");
		strcpy(dest->right_soft_clipped_qual , "");
		dest->left_soft_clipped_sequence_length = 0;
		dest->right_soft_clipped_sequence_length = 0;
		strcpy(dest->soft_clips_removed_sequence , dest->sequence);
		strcpy(dest->soft_clips_removed_quality_scores , dest->quality_scores);
		dest->soft_clips_removed_sequence_len = 0;
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

unsigned short int readSingleAlignmentFromSAMAlignmentFile (
		FILE *fhr,
		struct Sam_Alignment *s,
		char *ended,
		char *alignment_format,
		unsigned short int AS_tag_presence,
		char **split_line,
		char **split_tags)
{
	/*************************************************************************************************************************
	 * Reads in each alignment and stores the values in them Sam_Alignment object
	 * Decided which fields to retain based on whether the alignment was SE or PE
	 **************************************************************************************************************************/

	size_t len = 0;
	ssize_t line_len;

	char *line = NULL; // for reading each line
	unsigned short int number_of_fields;

	if ( strcmp (alignment_format , "SAM") == 0 )
	{
		line_len = getline ( &line , &len , fhr);
		if ( line_len == -1 ) return -1;
		number_of_fields = splitByDelimiter (line , '\t' , split_line);
		populateSamAlignmentInstance (s ,
				split_line ,
				split_tags ,
				number_of_fields ,
				AS_tag_presence);
	}
	else if ( strcmp (alignment_format , "BAM") == 0 )
	{

	}
	return 0;
}

#endif /* ABRIDGE_FUNCTIONS_DEFINITIONS_H_ */
