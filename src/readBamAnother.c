/* prints bases and quals */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <bam.h>
#define WHERE fprintf(stderr, "%d\n", __LINE__)
int main(int argc, char **argv)
{
    bam1_t *b = bam_init1();
    bamFile in = bam_open(argv[1], "r");
    bam_header_t *header;
    if (in == NULL)
        return -1;
    if (b == NULL)
        return -1;
    header = bam_header_read(in);
    while (bam_read1(in, b) >= 0)
    {
        int i;
        const bam1_core_t *c = &b->core;
        uint8_t *s = bam1_seq(b), *t = bam1_qual(b);

        fwrite(bam1_qname(b), c->l_qname - 1, sizeof(char), stdout);
        fputc('\t', stdout);
        for (i = 0; i < c->l_qseq; ++i)
            fputc(bam_nt16_rev_table[bam1_seqi(s, i)], stdout);
        fputc('\t', stdout);
        if (t[0] == 0xff)
        {
            fputs("*", stdout);
        }
        else
        {
            for (i = 0; i < c->l_qseq; ++i)
                fputc(t[i] + 33, stdout);
        }
        fputc('\n', stdout);
    }
    bam_header_destroy(header);
    bam_close(in);
    bam_destroy1(b);
    return 0;
}