# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <math.h>
# include <stdbool.h>
# include "data_structure_definitions.h"
# include "function_definitions.h"

void generateCoverageFromCompressedMappedFile (
		char *pass1_filename,
		char *abridge_index_filename,
		int d,
		int bg,
		int bga,
		int split,
		int generate_overlapping_coverage,
		int generate_nonoverlapping_coverage,
		int single,
		int max_reads_in_each_line,
		char *dictionary_name)
{
	/****************************************************************************************************************************************
	 * Variable declaration
	 ****************************************************************************************************************************************/
	char chromosome[1000];
	char samformatflag_replacer_characters[] = "BEFHIJKLMOPQRSUVWXYZbdefhijklmopqrsuvwxyz";

	char *temp;
	char *last_cigar;
	char *convert_to_int_temp;
	char *buffer_for_pass1;
	char *buffer_for_index;
	char *ptr_to_icigars;
	char **split_on_newline;
	char **split_on_tab;
	char **split_on_dash;
	char **split_on_comma;
	char **split_on_delimiter;
	char **split_on_newline_qual;
	char **read_names;

	size_t len = 0;
	ssize_t line_len;

	int number_of_bytes_read_from_compressed_file;
	int max_bytes_for_current_index_entry;
	int i, j, k, l;
	int number_of_entries_in_cluster;
	int number_of_elements_after_split_on_delimiter;
	int BUFFER_SIZE = 8 * 100 * 1024 * 1024; // 100 MB
	//int ROWS_split_on_newline = ROWS * 10; //10,000
	//int COLS_split_on_newline = COLS * 1000; //1,000,000
	int ROWS_split_on_tab = 10; //10
	int COLS_split_on_tab = COLS * 10; //100,000
	int ROWS_split_on_dash = 5; //5
	int COLS_split_on_dash = MAX_SEQ_LEN * 3; //3,000
	int ROWS_split_on_comma = ROWS * 10; //10,000
	int COLS_split_on_comma = MAX_SEQ_LEN * 3; //3,000
	int ROWS_split_on_newline = ROWS * 100; // 1,000,000
	int COLS_split_on_newline = COLS * 100; // 1,000,000
	int number_of_unique_samformatflags;
	int length_of_continuous_segment;
	int number_of_commas_in_each_line;
	int max_number_of_commas;
	int number_of_distinct_cigars_in_a_line;
	int number_of_repititions_of_the_same_reads;
	int num_of_cigar_types;

	long long int *coverage_array;
	long long int curr_position;
	long long int prev_stopping_location;
	long long int last_location_of_current_chromosome;
	long long int splice_distance;

	short int flag_ignore_mismatches;
	short int flag_ignore_soft_clippings;
	short int flag_ignore_unmapped_sequences;
	short int flag_ignore_quality_score;
	short int flag_save_all_quality_scores;
	short int flag_save_exact_quality_scores;
	short int flag_save_scores;
	short int number_of_columns;
	short int first_position;

	FILE *fhr_pass1;
	FILE *fhr_index;
	FILE *fhr_dictionary;

	struct Abridge_Index *abridge_index;
	struct Cigar_Items cigar_items_instance[MAX_SEQ_LEN];
	/****************************************************************************************************************************************/

	/****************************************************************************************************************************************
	 * Variable initialization
	 ****************************************************************************************************************************************/
	fhr_pass1 = fopen (pass1_filename , "rb");
	if ( fhr_pass1 == NULL )
	{
		printf ("Error! File pass1_filename = %s not found" , pass1_filename);
		exit (1);
	}
	fhr_index = fopen (abridge_index_filename , "r");
	if ( fhr_index == NULL )
	{
		printf ("Error! File abridge_index_filename = %s not found" , abridge_index_filename);
		exit (1);
	}
	if ( single == 0 )
	{
		fhr_dictionary = fopen (dictionary_name , "r");
		if ( fhr_dictionary == NULL )
		{
			printf ("Error! File dictionary_name = %s not found" , dictionary_name);
			exit (1);
		}
	}

	split_on_tab = ( char** ) malloc (sizeof(char*) * ROWS_split_on_tab);
	for ( i = 0 ; i < ROWS_split_on_tab ; i++ )
		split_on_tab[i] = ( char* ) malloc (sizeof(char) * COLS_split_on_tab);

	max_reads_in_each_line += 10;
	read_names = ( char** ) malloc (sizeof(char*) * max_reads_in_each_line);
	for ( i = 0 ; i < max_reads_in_each_line ; i++ )
		read_names[i] = ( char* ) malloc (sizeof(char) * 100);

	split_on_dash = ( char** ) malloc (sizeof(char*) * ROWS_split_on_dash);
	for ( i = 0 ; i < ROWS_split_on_dash ; i++ )
		split_on_dash[i] = ( char* ) malloc (sizeof(char) * COLS_split_on_dash);

	split_on_comma = ( char** ) malloc (sizeof(char*) * ROWS_split_on_comma);
	for ( i = 0 ; i < ROWS_split_on_comma ; i++ )
		split_on_comma[i] = ( char* ) malloc (sizeof(char) * COLS_split_on_comma);

	split_on_newline = ( char** ) malloc (sizeof(char*) * ROWS_split_on_newline);
	for ( i = 0 ; i < ROWS_split_on_comma ; i++ )
		split_on_newline[i] = ( char* ) malloc (sizeof(char) * COLS_split_on_newline);

	abridge_index = allocateMemoryAbridge_Index ();
	last_cigar = ( char* ) malloc (sizeof(char) * 1000);
	last_cigar[0] = '\0';
	if ( strcmp (dictionary_name , "dummy") == 0 ) max_reads_in_each_line = -1;
	/****************************************************************************************************************************************/

	line_len = getline ( &buffer_for_index , &len , fhr_index);
	splitByDelimiter (buffer_for_index , '\t' , split_on_tab);

	flag_ignore_mismatches = strtol (split_on_tab[0] , &convert_to_int_temp , 10);
	flag_ignore_soft_clippings = strtol (split_on_tab[1] , &convert_to_int_temp , 10);
	flag_ignore_unmapped_sequences = strtol (split_on_tab[2] , &convert_to_int_temp , 10);
	flag_ignore_quality_score = strtol (split_on_tab[3] , &convert_to_int_temp , 10);
	flag_save_all_quality_scores = strtol (split_on_tab[4] , &convert_to_int_temp , 10);
	flag_save_exact_quality_scores = strtol (split_on_tab[5] , &convert_to_int_temp , 10);
	flag_save_scores = strtol (split_on_tab[6] , &convert_to_int_temp , 10);

	readAbridgeIndex (abridge_index , abridge_index_filename , split_on_newline , &flag_ignore_mismatches , &flag_ignore_soft_clippings , &flag_ignore_unmapped_sequences , &flag_ignore_quality_score , &flag_save_all_quality_scores , &flag_save_exact_quality_scores , &flag_save_scores);

	line_len = getline ( &buffer_for_pass1 , &len , fhr_pass1); // Reading the first line and getting rid of it.
	coverage_array = NULL;
	for ( i = 0 ; i < abridge_index->number_of_items ; i++ )
	{
		length_of_continuous_segment = abridge_index->end[i] - abridge_index->start[i] + 1;
		max_bytes_for_current_index_entry = abridge_index->end_byte[i] - abridge_index->start_byte[i];
		number_of_bytes_read_from_compressed_file = 0;
		if ( coverage_array != NULL ) free (coverage_array);
		coverage_array = ( long long int* ) malloc (sizeof(long long int*) * length_of_continuous_segment);
		for ( j = 0 ; j < length_of_continuous_segment ; j++ )
			coverage_array[j] = 0;
		while ( number_of_bytes_read_from_compressed_file < max_bytes_for_current_index_entry )
		{
			/*
			 printf ("\n%s number_of_bytes_read_from_compressed_file = %d max_bytes_for_current_index_entry = %d" , abridge_index->chromosome[i] , number_of_bytes_read_from_compressed_file , max_bytes_for_current_index_entry);
			 fflush (stdout);
			 */
			line_len = getline ( &buffer_for_pass1 , &len , fhr_pass1);
			if ( buffer_for_pass1[0] == '@' )
			{
				/***********************************************************************************
				 * Printing out the last zero coverage nucleotides for the previous chromosome
				 * The final nucleotides of the last chromosome has been printed outside the loop
				 ***********************************************************************************/
				if ( d == 1 && bg == 0 && bga == 0 )
				{
					if ( i != 0 ) // Not the first chromosome
					{
						//printf ("\nThis %s start = %d end = %d" , abridge_index->chromosome[i - 1] , prev_stopping_location + 1 , last_location_of_current_chromosome);
						for ( j = prev_stopping_location + 1 ;
								j <= last_location_of_current_chromosome ; j++ )
							printf ("\n%s\t%d\t%d" , abridge_index->chromosome[i - 1] , j , 0);
					}
				}
				else if ( d == 0 && bg == 1 && bga == 0 )
				{
					/*
					 * Nucleotides with 0 coverage is skipped
					 */
				}
				else if ( d == 0 && bg == 0 && bga == 1 )
				{
					if ( i != 0 && prev_stopping_location + 1 - 1 != last_location_of_current_chromosome ) // Not the first chromosome
						printf ("\n%s\t%d\t%d\t%d" , abridge_index->chromosome[i - 1] , prev_stopping_location + 1 - 1 , last_location_of_current_chromosome , 0);
				}
				/***********************************************************************************/

				curr_position = 0;
				first_position = 1;
				prev_stopping_location = 0;
				number_of_columns = splitByDelimiter (buffer_for_pass1 , '\t' , split_on_tab);
				splitByDelimiter (split_on_tab[2] , ':' , split_on_comma);
				last_location_of_current_chromosome = strtol (split_on_comma[1] , &temp , 10);
				continue;
			}
			number_of_bytes_read_from_compressed_file += line_len;
			number_of_commas_in_each_line = 0;
			for ( j = 0 ; buffer_for_pass1[j] != '\0' ; j++ )
				if ( buffer_for_pass1[j] == ',' )
					number_of_commas_in_each_line++;
			if ( max_number_of_commas < number_of_commas_in_each_line )
				max_number_of_commas = number_of_commas_in_each_line;

			//printf ("\nline_len %d len %d" , line_len , len);
			if ( line_len > COLS_split_on_tab )
			{
				//printf ("\nB--> line_len %d COLS_split_on_tab %d" , line_len , COLS_split_on_tab);
				//fflush (stdout);
				for ( j = 0 ; j < ROWS_split_on_tab ; j++ )
					free (split_on_tab[j]);
				COLS_split_on_tab = line_len + 100;
				for ( j = 0 ; j < ROWS_split_on_tab ; j++ )
					split_on_tab[j] = ( char* ) malloc (sizeof(char) * COLS_split_on_tab);
				//printf ("\nA--> line_len %d COLS_split_on_tab %d" , line_len , COLS_split_on_tab);
				//fflush (stdout);

			}
			if ( max_number_of_commas > ROWS_split_on_comma )
			{
				//printf ("\nB--> max_number_of_commas %d ROWS_split_on_comma %d" , max_number_of_commas , ROWS_split_on_comma);
				//fflush (stdout);
				for ( j = 0 ; j < ROWS_split_on_comma ; j++ )
					free (split_on_comma[i]);
				free (split_on_comma);
				ROWS_split_on_comma = line_len / max_number_of_commas + 10;
				split_on_comma = ( char** ) malloc (sizeof(char*) * ROWS_split_on_comma);
				for ( j = 0 ; j < ROWS_split_on_comma ; j++ )
					split_on_comma[j] = ( char* ) malloc (sizeof(char) * COLS_split_on_comma);
				//printf ("\nA--> max_number_of_commas %d ROWS_split_on_comma %d" , max_number_of_commas , ROWS_split_on_comma);
				//fflush (stdout);
			}
			number_of_columns = splitByDelimiter (buffer_for_pass1 , '\t' , split_on_tab);
			if ( max_reads_in_each_line == -1 ) number_of_columns++;
			//printf ("\nmax_reads_in_each_line %d dictionary name %s" , max_reads_in_each_line , dictionary_name);
			//fflush (stdout);
			if ( number_of_columns == 2 )
			{
				curr_position++;
				ptr_to_icigars = split_on_tab[0];
			}
			else if ( number_of_columns == 3 )
			{
				curr_position += strtol (split_on_tab[0] , &convert_to_int_temp , 10);
				ptr_to_icigars = split_on_tab[1];
			}
			//printf ("\nChromosome=%s curr_position=%d" , abridge_index->chromosome[i] , curr_position);
			/*
			 if ( abridge_index->start[i] != 23162 || abridge_index->end[i] != 207509 )
			 continue;
			 */
			number_of_distinct_cigars_in_a_line = splitByDelimiter (ptr_to_icigars , ',' , split_on_comma);

			for ( j = 0 ; j < number_of_distinct_cigars_in_a_line ; j++ )
			{
				splitByDelimiter (split_on_comma[j] , '-' , split_on_dash);
				number_of_repititions_of_the_same_reads = strtol (split_on_dash[1] , &temp , 10);
				if ( split_on_comma[j][1] == '-' && isalpha (split_on_dash[0][0]) != 0 )
				{
					//printf ("\ncigar Entering here borrowing from %s" , last_cigar);
					//printf ("\ncigar %s %d" , ptr_to_icigars , number_of_distinct_cigars_in_a_line);
					splitByDelimiter (last_cigar , '-' , split_on_dash);
					splitCigar (split_on_dash[0] , &num_of_cigar_types , cigar_items_instance);
				}
				else
				{
					splitCigar (split_on_dash[0] , &num_of_cigar_types , cigar_items_instance);
					strcpy(last_cigar , split_on_dash[0]);
				}
				//printf ("\nnum_of_cigar_types %d" , num_of_cigar_types);
				l = 0;
				//splice_distance = 0;
				for ( k = 0 ; k < num_of_cigar_types ; k++ )
				{

					if ( cigar_items_instance[k].def == 'a' || cigar_items_instance[k].def == 't' || cigar_items_instance[k].def == 'g' || cigar_items_instance[k].def == 'c' || cigar_items_instance[k].def == 'n' ) // Soft clips
						continue;
					else if ( cigar_items_instance[k].def >= ( 33 + 90 ) && cigar_items_instance[k].def <= ( 73 + 90 ) ) // Quality scores
						continue;
					else if ( cigar_items_instance[k].def == 'N' && split == 1 ) // Intron splice
					{
						l += cigar_items_instance[k].len;
						continue;
					}
					else if ( cigar_items_instance[k].def == 'D' )
					{
						l += cigar_items_instance[k].len;
						continue;
					}
					else if ( cigar_items_instance[k].def == '!' || cigar_items_instance[k].def == '"' || cigar_items_instance[k].def == '#' || cigar_items_instance[k].def == '$' || cigar_items_instance[k].def == '%' ) // Insertions in reads
						continue;
					else if ( cigar_items_instance[k].def == '&' || cigar_items_instance[k].def == '\'' || cigar_items_instance[k].def == '(' || cigar_items_instance[k].def == ')' || cigar_items_instance[k].def == '*' )
						cigar_items_instance[k].len = 1;
					/*
					 if ( abridge_index->start[i] == 23046848 )
					 printf ("\ncigar %d %c %d" , cigar_items_instance[k].len , cigar_items_instance[k].def , curr_position);
					 */
					if ( generate_overlapping_coverage == 1 && generate_nonoverlapping_coverage == 0 )
					{
						//printf ("\n%d %c" , cigar_items_instance[k].len , cigar_items_instance[k].def);
						while ( cigar_items_instance[k].len-- )
						{
							//printf ("\nUpdating coverage");
							if ( curr_position - abridge_index->start[i] + l < 0 || curr_position - abridge_index->start[i] + l >= length_of_continuous_segment )
							{
								printf ("\nIncorrect");
								printf ("\n%d %d" , curr_position - abridge_index->start[i] + l , length_of_continuous_segment);
								printf ("\nCurrent position=%d, abrdige_start=%d l=%d" , curr_position , abridge_index->start[i] , l);
								printf ("\niCIGAR: %s j = %d k = %d" , split_on_comma[j] , j , k);
								fflush (stdout);
								exit (1);
							}
							coverage_array[curr_position - abridge_index->start[i] + l] += number_of_repititions_of_the_same_reads;
							l++;
						}
					}
					else if ( generate_overlapping_coverage == 0 && generate_nonoverlapping_coverage == 1 )
						coverage_array[curr_position - abridge_index->start[i]] += number_of_repititions_of_the_same_reads;
				}
				/*
				 if ( abridge_index->start[i] == 23046848 )
				 printf ("\ncigar====================================================================================================================");
				 */
			}
		}
		/*
		 printf ("\n%s number_of_bytes_read_from_compressed_file = %d max_bytes_for_current_index_entry = %d" , abridge_index->chromosome[i] , number_of_bytes_read_from_compressed_file , max_bytes_for_current_index_entry);
		 printf ("\n=======================================================================================================================================");
		 fflush (stdout);
		 */
		/*
		 * Print the coverage as requested by user
		 */
		if ( d == 1 && bg == 0 && bga == 0 )
		{
			if ( first_position == 1 )
			{
				for ( j = 1 ; j < abridge_index->start[i] ; j++ )
					printf ("\n%s\t%d\t%d" , abridge_index->chromosome[i] , j , 0);
				first_position = 0;
			}
			else if ( first_position == 0 )
			{
				for ( j = prev_stopping_location + 1 ;
						j < abridge_index->start[i] ; j++ )
					printf ("\n%s\t%d\t%d" , abridge_index->chromosome[i] , j , 0);

			}
			for ( j = 0 ; j < length_of_continuous_segment ; j++ )
				printf ("\n%s\t%d\t%d" , abridge_index->chromosome[i] , abridge_index->start[i] + j , coverage_array[j]);
			prev_stopping_location = abridge_index->start[i] + j - 1;
		}
		else if ( d == 0 && ( bg == 1 || bga == 1 ) )
		{
			if ( first_position == 1 )
			{
				if ( 1 - 1 != abridge_index->start[i] - 1 && bg == 0 && bga == 1 )
					printf ("\n%s\t%d\t%d\t%d" , abridge_index->chromosome[i] , 1 - 1 , abridge_index->start[i] - 1 , 0);
				first_position = 0;
			}
			else if ( first_position == 0 && bg == 0 && bga == 1 )
			{
				if ( prev_stopping_location + 1 - 1 != abridge_index->start[i] - 1 )
					printf ("\n%s\t%d\t%d\t%d" , abridge_index->chromosome[i] , prev_stopping_location + 1 - 1 , abridge_index->start[i] - 1 , 0);
			}
			int local_start, local_end, previous_coverage;
			int print_last_record = 0;
			for ( j = 0 ; j < length_of_continuous_segment ; j++ )
			{
				if ( j == 0 )
				{
					previous_coverage = coverage_array[j];
					local_start = 0;
					local_end = 0;
				}
				else
				{
					if ( previous_coverage == coverage_array[j] )
					{
						local_end++;
						if ( j == length_of_continuous_segment - 1 )
							print_last_record = 1;
					}
					else //Coverage mismatch
					{
						if ( abridge_index->start[i] + local_start - 1 != abridge_index->start[i] + local_end )
						{
							if ( ( bg == 1 && previous_coverage != 0 ) || bga == 1 )
								printf ("\n%s\t%d\t%d\t%d" , abridge_index->chromosome[i] , abridge_index->start[i] + local_start - 1 , abridge_index->start[i] + local_end , previous_coverage);
						}
						/*
						 if ( j == length_of_continuous_segment - 1 )
						 printf ("Last nucleotide group of cluster");
						 */
						previous_coverage = coverage_array[j];
						if ( j != length_of_continuous_segment - 1 )
						{
							local_start = local_end + 1;
							local_end++;
						}
						else print_last_record = 0;
					}
				}
				//printf ("\n%s\t%d\t%d" , abridge_index->chromosome[i] , abridge_index->start[i] + j , coverage_array[j]);
			}
			if ( print_last_record == 1 )
			{
				if ( abridge_index->start[i] + local_start - 1 != abridge_index->start[i] + local_end )
					if ( ( bg == 1 && previous_coverage != 0 ) || bga == 1 )
						printf ("\n%s\t%d\t%d\t%d" , abridge_index->chromosome[i] , abridge_index->start[i] + local_start - 1 , abridge_index->start[i] + local_end , previous_coverage);
			}
			if ( coverage_array[length_of_continuous_segment - 1] != coverage_array[length_of_continuous_segment - 2] ) // The last nucleotide was a mismatch
			{
				if ( abridge_index->start[i] + length_of_continuous_segment - 2 != abridge_index->start[i] + length_of_continuous_segment - 1 )
					if ( ( bg == 1 && previous_coverage != 0 ) || bga == 1 )
						printf ("\n%s\t%d\t%d\t%d" , abridge_index->chromosome[i] , abridge_index->start[i] + length_of_continuous_segment - 2 , abridge_index->start[i] + length_of_continuous_segment - 1 , previous_coverage);
			}
			prev_stopping_location = abridge_index->start[i] + j - 1;
		}
		//printf (" Cluster end Previous stopping location %d cluster start %d length of continuous segment %d" , prev_stopping_location , abridge_index->start[i] , length_of_continuous_segment);
	}

	/*
	 * Print out the zero entries for the final nucleotides of the last chromosome
	 */
	if ( d == 1 && bg == 0 && bga == 0 )
	{
		for ( j = prev_stopping_location + 1 ;
				j <= last_location_of_current_chromosome ; j++ )
			printf ("\n%s\t%d\t%d" , abridge_index->chromosome[i] , j , 0);
	}
	else if ( d == 0 && bg == 1 && bga == 0 )
	{
		/*
		 * Nucleotides with 0 coverage is skipped
		 */
	}
	else if ( d == 0 && bg == 0 && bga == 1 )
	{
		if ( prev_stopping_location + 1 - 1 != last_location_of_current_chromosome )
			printf ("\n%s\t%d\t%d\t%d" , abridge_index->chromosome[i - 1] , prev_stopping_location + 1 - 1 , last_location_of_current_chromosome , 0);
	}
	fflush (stdout);
}

int main (int argc, char *argv[])
{
	/****************************************************************************************************************************************
	 * Variable declaration
	 ****************************************************************************************************************************************/
	char pass1_filename[FILENAME_LENGTH];
	char abridge_index_filename[FILENAME_LENGTH];
	char dictionary_name[FILENAME_LENGTH];

	char *temp;

	int d;
	int bg;
	int bga;
	int split;
	int generate_overlapping_coverage;
	int generate_nonoverlapping_coverage;
	int single;
	int max_reads_in_each_line;
	/****************************************************************************************************************************************/

	/****************************************************************************************************************************************
	 * Variable initialization
	 ****************************************************************************************************************************************/
	strcpy(pass1_filename , argv[1]);
	strcpy(abridge_index_filename , argv[2]);
	d = strtol (argv[3] , &temp , 10);
	bg = strtol (argv[4] , &temp , 10);
	bga = strtol (argv[5] , &temp , 10);
	split = strtol (argv[6] , &temp , 10);
	generate_overlapping_coverage = strtol (argv[7] , &temp , 10);
	generate_nonoverlapping_coverage = strtol (argv[8] , &temp , 10);
	single = strtol (argv[9] , &temp , 10);
	if ( argc > 10 )
	{
		strcpy(dictionary_name , argv[10]);
		max_reads_in_each_line = strtol (argv[11] , &temp , 10);
	}
	else
	{
		strcpy(dictionary_name , "dummy");
		max_reads_in_each_line = -1;
	}
	/****************************************************************************************************************************************/

	generateCoverageFromCompressedMappedFile (pass1_filename , abridge_index_filename , d , bg , bga , split , generate_overlapping_coverage , generate_nonoverlapping_coverage , single , max_reads_in_each_line , dictionary_name);

	return 0;
}
