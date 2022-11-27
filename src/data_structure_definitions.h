#ifndef _ABRIDGE_DATA_STRUCTURE_DEFINITIONS_H_
#define _ABRIDGE_DATA_STRUCTURE_DEFINITIONS_H_

/*************************************************************************************************************************
 * Define macros
 **************************************************************************************************************************/
#define MAX_FILENAME_LENGTH 10000
#define TEN 10
#define ONE_HUNDRED 100
#define ONE_THOUSAND 1000
#define TEN_THOUSAND 10000
#define HUNDRED_THOUSAND 100000
#define MAX_SEQ_LEN 1000
#define MAX_CIGAR_ITEMS 5000
#define MAX_CIGAR_LENGTH 10000
#define MAX_TAG_LENGTH 1000
#define MAX_POOL_SIZE 1000000
#define MIN_POOL_SIZE 10000
#define MAX_UNIQUE_CIGARS 10000
#define MAX_LINE_LENGTH_TO_BE_WRITTEN_TO_FILE 10000000
#define MAX_LENGTH_RLE 1000
#define MAX_CIGAR_FREQ_SIZE 1000000
#define MAX_ICIGAR_LENGTH 100000
#define MAX_ICIGAR_LENGTH_PASS1_COL2 10000000
#define MAX_REFERENCE_SEQUENCES 300000
#define MAX_GENERAL_LEN 1000000
#define MAX_BUFFER_SIZE_FOR_READING_PASS2_FILE 1073741824
#define MAX_REFERENCE_SEQ_LEN 1000000000
#define MAX_FILES_FOR_MERGING 1000
#define MAX_READ_ID_LENGTH 1000
#define QUAL_SCORE_ADJUSTMENT 100

/*************************************************************************************************************************
 * Define characters for mismatch and deletion
 *************************************************************************************************************************/

// ATGCN - ASCII Codes 33-37
char insert_characters[] =
    {'!', '"', '#', '$', '%', '\0'};
// ATGCN - ASCII Codes 38-42
char mismatch_characters[] =
    {'&', '\'', '(', ')', '*', '\0'};

struct Reference_Sequence_Info
{
    /*************************************************************************************************************************
     * Data structure definition to outline the how nucleotides from the reference will be loaded into the memory
     **************************************************************************************************************************/
    char *reference_name;                   // Name of the reference, derived from the fasta file
    char *nucleotide_sequence;              // Should contain only ATGCN
    unsigned long long int sequence_length; // Length of the sequence - should be calculated only once
};

struct Reference
{
    /*************************************************************************************************************************
     * Data structure definition to outline how the reference information will be organized in the memory
     **************************************************************************************************************************/
    unsigned int number_of_reference_sequences;     // Total number of reference sequences - to be inferred from the file
    struct Reference_Sequence_Info **each_sequence; // Data structure for each sequence
};

struct Cigar_Items
{
    /*************************************************************************************************************************
     * Data structure definition to outline how CIGAR items are loaded and stored
     *************************************************************************************************************************/
    char def;
    int len;
};

struct Soft_Clippings
{
    /*************************************************************************************************************************
     * Data structure definition to outline how soft clips from each alignment is stored
     **************************************************************************************************************************/
    char *left_sequence;                // Nucleotide sequence soft clipped from the left
    char *right_sequence;               // Nucleotide sequence soft clipped from the right
    char *left_qual;                    // Quality scores soft clipped from the left
    char *right_qual;                   // Quality scores soft clipped from the right
    unsigned int left_sequence_length;  // Store the length of the left clipped sequence
    unsigned int right_sequence_length; // Store the length of the right clipped sequence
    char *soft_clips_removed_seq;       // Nucleotides from the portion of the sequence not soft clipped
    char *soft_clips_removed_qual;      // Quality scores from the portion of the sequence not soft clipped
    int soft_clips_removed_seq_len;     // Length of the sequence without any soft clips
};

struct Sam_Alignment
{
    /*************************************************************************************************************************
     * Data structure definition to outline how each alignment will be loaded in the memory
     **************************************************************************************************************************/
    // Fields are listed in the order of appearance in SAM file
    /*
     * Mandatory fields
     */
    char *read_name;                       // Name of the read
    char *samflag;                         // Flag produced by aligner. For more details see https://www.samformat.info/sam-format-flag
    char *reference_name;                  // Chromosome name
    unsigned long long int start_position; // starting position of the read
    char *mapping_quality_score;           // mapping score - not very important but might be of use to some downstream software
    char *cigar;                           // the most important bit of information
    char *reference_name_next_mate;        // Reference name where the mate/next read is mapped
    char *start_position_next;             // Position of the mate/next read
    char *template_length;                 // Observed template length
    char *sequence;                        // The nucleotide read sequence
    char *quality_scores;                  // The read quality scores

    /*
     * Tags important for ABRIDGE compression
     */
    char *NH; // Stores the NH value, will be needed for generating iCIGAR
    char *AS; // Stores the AS value, will be needed for generating iCIGAR
    char *MD; // Stores the MD value, will be needed for generating iCIGAR

    /*
     * Other fields created by abridge
     */
    int read_seq_len;                     // basically strlen(seq). Field is kept to prevent extra calls to strlen
    struct Soft_Clippings soft_clippings; // Stores the left and the right soft clippings
    struct Cigar_Items *cigar_items;
    int number_of_cigar_items;

    char *cigar_extended; // Required for generating iCIGAR
    char *md_extended;    // Required for generating iCIGAR
    char *icigar;         // Stores the integrated representation comprising of all relevant information about the alignment - the iCIGAR
    char **splices;
};

struct Paired_Ended_Dictionary_Items
{
    /*************************************************************************************************************************
     * Data structure definition to outline how the dictionary is constructed for PE reads
     **************************************************************************************************************************/
    char *character;
    int *samflags;
};

#endif /* _ABRIDGE_DATA_STRUCTURE_DEFINITIONS_H_ */