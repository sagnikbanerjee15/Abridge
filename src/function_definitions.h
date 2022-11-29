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

struct Reference_Sequence_Info *allocateMemoryReference_Sequence_Info()
{
    /*************************************************************************************************************************
     * Data structure memory allocation function for Reference_Sequence_Info
     * Initialized memory only for the data structure itself and not for the underlying data
     * Memory will be dynamically allocated for the data structure variables during file read operation
     **************************************************************************************************************************/
    struct Reference_Sequence_Info *s;
    s = (struct Reference_Sequence_Info *)malloc(sizeof(struct Reference_Sequence_Info));
    s->reference_name = NULL;
    s->nucleotide_sequence = NULL;
    return s;
}

struct Reference *allocateMemoryReference(unsigned int n)
{
    /*************************************************************************************************************************
     * Data structure memory allocation function for Reference
     *************************************************************************************************************************/
    struct Reference *s;
    unsigned int i;

    s = (struct Reference *)malloc(sizeof(struct Reference));
    s->number_of_reference_sequences = n;
    s->each_sequence = (struct Reference_Sequence_Info **)malloc(
        sizeof(struct Reference_Sequence_Info *) * s->number_of_reference_sequences);
    for (i = 0; i < s->number_of_reference_sequences; i++)
    {
        s->each_sequence[i] = allocateMemoryReference_Sequence_Info();
    }
    return s;
}

struct Sam_Alignment *allocateMemorySam_Alignment()
{
    /********************************************************************
     * Allocates memory for SAM alignments
     ********************************************************************/
    /********************************************************************
     * Variable declarations
     ********************************************************************/
    struct Sam_Alignment *s;
    int i;
    /********************************************************************/
    s = (struct Sam_Alignment *)malloc(sizeof(struct Sam_Alignment));

    s->read_name = (char *)malloc(sizeof(char) * MAX_SEQ_LEN);
    s->samflag = (char *)malloc(sizeof(char) * TEN);
    s->reference_name = (char *)malloc(sizeof(char) * MAX_SEQ_LEN);
    s->start_position = 0;
    s->mapping_quality_score = (char *)malloc(sizeof(char) * TEN);
    s->cigar = (char *)malloc(sizeof(char) * MAX_SEQ_LEN);
    s->reference_name_next_mate = (char *)malloc(sizeof(char) * MAX_SEQ_LEN);
    s->start_position_next = (char *)malloc(sizeof(char) * MAX_SEQ_LEN);
    s->template_length = (char *)malloc(sizeof(char) * MAX_SEQ_LEN);
    s->sequence = (char *)malloc(sizeof(char) * MAX_SEQ_LEN);
    s->quality_scores = (char *)malloc(sizeof(char) * MAX_SEQ_LEN);

    s->soft_clippings.left_sequence = (char *)malloc(sizeof(char) * MAX_SEQ_LEN);
    s->soft_clippings.right_sequence = (char *)malloc(sizeof(char) * MAX_SEQ_LEN);
    s->soft_clippings.left_qual = (char *)malloc(sizeof(char) * MAX_SEQ_LEN);
    s->soft_clippings.right_qual = (char *)malloc(sizeof(char) * MAX_SEQ_LEN);
    s->soft_clippings.left_sequence_length = 0;
    s->soft_clippings.right_sequence_length = 0;
    s->soft_clippings.soft_clips_removed_qual = (char *)malloc(sizeof(char) * MAX_SEQ_LEN);
    s->soft_clippings.soft_clips_removed_seq = (char *)malloc(sizeof(char) * MAX_SEQ_LEN);
    s->soft_clippings.soft_clips_removed_seq_len = 0;

    s->cigar_items = (struct Cigar_Items *)malloc(sizeof(struct Cigar_Items) * MAX_CIGAR_ITEMS);
    s->number_of_cigar_items = 0;

    s->cigar_extended = (char *)malloc(sizeof(char) * (MAX_SEQ_LEN * 2));
    s->md_extended = (char *)malloc(sizeof(char) * (MAX_SEQ_LEN * 2));
    s->icigar = (char *)malloc(sizeof(char) * (MAX_SEQ_LEN * 2));
    s->splices = (char **)malloc(sizeof(char *) * 100);
    for (i = 0; i < 100; i++)
        s->splices[i] = (char *)malloc(sizeof(char) * 50);
    s->cigar_extended[0] = '\0';
    s->icigar[0] = '\0';
    s->md_extended[0] = '\0';

    s->NH = (char *)malloc(sizeof(char) * TEN);
    s->AS = (char *)malloc(sizeof(char) * TEN);
    s->MD = (char *)malloc(sizeof(char) * (MAX_SEQ_LEN * 2));

    strcpy(s->NH, "-1");
    strcpy(s->AS, "X");
    strcpy(s->MD, "-1");
    return s;
}

/*************************************************************************************************************************
 *************************************************************************************************************************
 * Functions for allocating memory to Data Structures defined in data_structure_definitions.h
 *************************************************************************************************************************
 * END
 *************************************************************************************************************************/

char *strupr(char *s)
{
    /*************************************************************************************************************************
     * Convert string to uppercase
     *************************************************************************************************************************/
    char *pc = s;
    while (*pc)
    {
        *pc = toupper(*pc); /* toupper() requires <ctype.h> */
        ++pc;
    }
    return s;
}

void readCompleteReference(
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

    fhr = fopen(reference_filename, "rb");
    if (fhr == NULL)
    {
        printf("Error! File not found %s", reference_filename);
        exit(1);
    }

    /*************************************************************************************************************************
     * 1st pass to determine number of reference sequences
     *************************************************************************************************************************/
    while ((line_len = getline(&buffer, &len, fhr)) != -1)
        if (buffer[0] == '>')
            number_of_reference_sequences++;

    /*************************************************************************************************************************
     * Allocate memory for the reference
     *************************************************************************************************************************/
    reference = allocateMemoryReference(number_of_reference_sequences);

    /*************************************************************************************************************************
     * 2nd pass for loading all information from the fasta file
     *************************************************************************************************************************/
    rewind(fhr);
    number_of_reference_sequences_index = -1;

    while ((line_len = getline(&buffer, &len, fhr)) != -1)
    {
        // printf("\n%lld", strlen(buffer));
        if (strlen(buffer) <= 1)
            continue;
        if (buffer[0] == '>')
        {
            number_of_reference_sequences_index++;
            reference->each_sequence[number_of_reference_sequences_index]->reference_name = (char *)malloc(sizeof(char) * (line_len + 1));
            /*************************************************************************************************************************
             * Copy the name of the reference only till the first space
             *************************************************************************************************************************/
            for (i = 1, j = 0; buffer[i] != 32; i++, j++)
                reference->each_sequence[number_of_reference_sequences_index]->reference_name[j] = buffer[j];
            reference->each_sequence[number_of_reference_sequences_index]->reference_name[j] = '\0';
        }
        else
        {
            reference->each_sequence[number_of_reference_sequences_index]->nucleotide_sequence = (char *)malloc(sizeof(char) * (line_len + 1));
            buffer = strupr(buffer);
            strcpy(
                reference->each_sequence[number_of_reference_sequences_index]->nucleotide_sequence,
                buffer);
            reference->each_sequence[number_of_reference_sequences_index]->sequence_length = strlen(buffer);
        }
    }

    fclose(fhr);
}

inline char *convertSignedIntegerToString(long long int x)
{
    /*************************************************************************************************************************
     * Converts signed integer to string
     **************************************************************************************************************************/
    char str[256];
    sprintf(str, "%lld", x);
}

inline char *convertUnsignedIntegerToString(unsigned long long int x)
{
    /*************************************************************************************************************************
     * Converts unsigned integers to string
    (*************************************************************************************************************************/
    char str[256];
    sprintf(str, "%llu", x);
}

inline unsigned long long int convertStringToUnsignedInteger(char *str)
{
    /*************************************************************************************************************************
     * Converts string to unsigned integer
     **************************************************************************************************************************/
    char convertStringToUnsignedInteger_value[10];
    char *convertStringToUnsignedInteger_eptr;

    return strtoull(convertStringToUnsignedInteger_value, &convertStringToUnsignedInteger_eptr, 10);
}

inline long long int convertStringToSignedInteger(char *str)
{
    /*************************************************************************************************************************
     * Converts string to signed integer
     **************************************************************************************************************************/
    char convertStringToSignedInteger_value[10];
    char *convertStringToSignedInteger_eptr;

    return strtoull(convertStringToSignedInteger_value, &convertStringToSignedInteger_eptr, 10);
}

int splitByDelimiter(char *line, char delimiter, char **new_string)
{
    /*********************************************************************
     * Splits a string on a requested character and
     * returns the number of splits
     *********************************************************************/

    /********************************************************************
     * Variable declarations
     ********************************************************************/
    int i = 0, j = 0, ctr = 0;
    /********************************************************************/

    for (i = 0; line[i] != '\0'; i++)
    {
        // printf ("\nsplitByDelimiter i=%i" , i);
        //  if space or NULL found, assign NULL into new_string[ctr]
        if (line[i] == delimiter)
        {
            new_string[ctr][j] = '\0';
            // printf ("\nsplitByDelimiter ctr=%i" , ctr);
            ctr++; // for next word
            j = 0; // for next word, init index to 0
        }
        else if (line[i] == '\n')
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

#endif /* ABRIDGE_FUNCTIONS_DEFINITIONS_H_ */