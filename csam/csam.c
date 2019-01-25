#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <limits.h>
#include <unistd.h>
#include <math.h>
#include <wchar.h>
#include <locale.h>
#include <time.h>
#include <inttypes.h>
#include <sys/ioctl.h>
#include <locale.h>
#include <getopt.h>
#include <assert.h>

#include <zlib.h>
#include <bam.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/khash.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>

#define VERSION "0.1.0"
#define CONTEXT_SIZE 1

#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define max(X, Y)  ((X) < (Y) ? (Y) : (X))
#define basename(str) (strrchr(str, '/') ? strrchr(str, '/') + 1 : str)

const char base[] = "ACGT";
int8_t seq_comp_table[16] = {0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15};

typedef struct
{
	char *bam, *ref, *pfx;
	int max_is, min_is;
	int min_mq, min_qs;
} arg_t;

struct option long_options[] =
{
	{"bam", required_argument, 0, 'b'},
	{"ref", required_argument, 0, 'r'},
	{"pfx", required_argument, 0, 'p'},
	{"max_is", required_argument, 0, 'M'},
	{"min_is", required_argument, 0, 'm'},
	{"min_mq", required_argument, 0, 'Q'},
	{"min_qs", required_argument, 0, 'q'},
	{"version", no_argument, 0, 'v'},
	{"help", no_argument, 0, 'h'},
	{0, 0, 0, 0}
};

void usage(char *str);
void ini_arg(arg_t *arg);
static inline void nt_to_upper(char *nt);
char *get_cigar(bam1_t *b);
char *get_read(const bam1_t *b);
uint64_t tot_reads(const char *fn);
void progressbar(uint64_t cur, uint64_t tot, int win);
int seq2idx(const char q, const char *r);
int qes2idx(const char q, const char *r);
void idx2seq(int idx, char *q, char *r);
int get_qual(const bam1_t *b, char **qual_out);
char *parse_rg(const bam_hdr_t *hdr, const char *tag);
void put_matrics(int m[][256], const bam_hdr_t *hdr, const char *prefix);

int main(int argc, char *argv[])
{
	if (argc == 1) usage(argv[0]);
	arg_t *arg = calloc(1, sizeof(arg_t));
	int c = 0, opt_idx = 0;
	uint64_t required_args = 0;
	while ((c = getopt_long(argc, argv,"b:r:p:M:m:Q:q:vh", long_options, &opt_idx)) != -1)
	{
		switch (c)
		{
			case 'b':
				arg->bam = optarg;
				required_args |= 1<<0;
				break;
			case 'r':
				arg->ref = optarg;
				required_args |= 1<<1;
				break;
			case 'p':
				arg->pfx = optarg;
				break;
			case 'M':
				arg->max_is = atoi(optarg);
				break;
			case 'm':
				arg->min_is = atoi(optarg);
				break;
			case 'Q':
				arg->min_mq = atoi(optarg);
				break;
			case 'q':
				arg->min_qs = atoi(optarg);
				break;
			case 'v':
				fputs(VERSION, stderr);
				fputc('\n', stderr);
				return(EXIT_SUCCESS);
			case 'h':
				usage(argv[0]);
				return EXIT_SUCCESS;
			case ':':
				fprintf(stderr, "Option -%c requires an argument\n", optopt);
				return EXIT_FAILURE;
				break;
			case '?':
				fprintf(stderr, "Option -%c is undefined\n", optopt);
				return EXIT_FAILURE;
				break;
		}
	}
	int required_args_set = __builtin_popcount(required_args & ((1<<2) - 1));
	if (required_args_set != 2)
	{
		fprintf(stderr, "[ERROR] %d required argument%s unspecified\n", 2 - required_args_set, (2 - required_args_set > 1) ? "s" : "");
		usage(argv[0]);
	}
	setlocale(LC_ALL, "");
	int ret = 0;
	int metrics[4][256] = {0};
	char q, r[4] = {0};
	// open bam and load ref index
	struct winsize w;
	ioctl(0, TIOCGWINSZ, &w);
	uint64_t n = 0, N = 0, tt = tot_reads(arg->bam);
	BGZF *in = bgzf_open(arg->bam, "r"); assert(in);
	faidx_t *fai = fai_load(arg->ref); assert(fai);
	bam_hdr_t *hdr = bam_hdr_read(in);
	bam1_t *b = bam_init1();
	while (bam_read1(in, b) >= 0)
	{
		++n;
		if (n * 100 / tt == N)
		{
			progressbar(n, tt, w.ws_col);
			++N;
		}
		if (!(b->core.flag & BAM_FPAIRED)) continue;
		if (b->core.flag & (BAM_FUNMAP | BAM_FDUP | BAM_FSECONDARY | BAM_FQCFAIL))
			continue;
		if (b->core.qual < arg->min_mq)
			continue;
		if (abs(b->core.isize) < arg->min_is || abs(b->core.isize) > arg->max_is)
			continue;
		char *qseq = get_read(b);
		char *qual = NULL;
		if (get_qual(b, &qual) < 0)
		{
			fprintf(stderr, "Error getting quality\n");
			exit(EXIT_FAILURE);
		}
		uint32_t *cigar = bam1_cigar(b);
		int end = bam_calend(&b->core, cigar);
		char *rseq = faidx_fetch_seq(fai, hdr->target_name[b->core.tid], b->core.pos - 1, end + 1, &ret);
		char *rseqp = rseq;
		while (*rseqp) nt_to_upper(rseqp++);
		int rn = 0, qn = 0, on = b->core.n_cigar;
		for(int i = 0; i < on; ++i)
		{
			int32_t op = bam_cigar_op(cigar[i]);
			int32_t ol = bam_cigar_oplen(cigar[i]);
			switch(op)
			{
				case BAM_CMATCH:
				case BAM_CDIFF:
				case BAM_CEQUAL:
					for (int j = 0; j < ol; ++j)
					{
						if ((int)qual[qn + j] - 33 < arg->min_qs)
							continue;
						strncpy(r, rseq + rn + j, 3);
						if (strchr(r, 'N')) continue;
						q = qseq[qn + j];
						if (q == 'N') continue;
						++metrics[(((b->core.flag & BAM_FREAD1) ? 1 : 0) << 1) | ((b->core.flag & BAM_FREVERSE) ? 1 : 0)][seq2idx(q, r)];
					}
					rn += ol;
					qn += ol;
					break;
				case BAM_CINS:
				case BAM_CSOFT_CLIP:
					qn += ol;
					break;
				case BAM_CDEL:
				case BAM_CREF_SKIP:
					rn += ol;
					break;
				case BAM_CPAD:
				case BAM_CBACK:
				case BAM_CHARD_CLIP:
					break;
				default:
					break;
			}
		}
		free(qseq);
		free(qual);
		free(rseq);
	}
	put_matrics(metrics, hdr, arg->pfx);
	bam_destroy1(b);
	bam_hdr_destroy(hdr);
	fai_destroy(fai);
	bgzf_close(in);
	free(arg);
	return 0;
}

/* CON 1 1; 0 0->0,3
 * PRO 1 0; 0 1->1,2
 */
void put_matrics(int m[][256], const bam_hdr_t *hdr, const char *prefix)
{
	char prea[] = "SAMPLE_ALIAS\tLIBRARY\tREF_BASE\tALT_BASE\tCONTEXT\tPRO_REF_BASES\tPRO_ALT_BASES\tCON_REF_BASES\tCON_ALT_BASES\tERROR_RATE\tQSCORE\n";
    char bait[] = "SAMPLE_ALIAS\tLIBRARY\tREF_BASE\tALT_BASE\tCONTEXT\tFWD_CXT_REF_BASES\tFWD_CXT_ALT_BASES\tREV_CXT_REF_BASES\tREV_CXT_ALT_BASES\tFWD_ERROR_RATE\tREV_ERROR_RATE\tERROR_RATE\tQSCORE\n";
	kstring_t ks_pa = {0, 0, 0}, ks_bb = {0, 0, 0};
	ksprintf(&ks_pa, "%s.pre_adapter_detail_metrics.txt", prefix ? prefix : "sample");
	ksprintf(&ks_bb, "%s.bait_bias_detail_metrics.txt", prefix ? prefix : "sample");
	FILE *fp_pa = fopen(ks_pa.s, "w");
	assert(fp_pa);
	FILE *fp_bb = fopen(ks_bb.s, "w");
	assert(fp_bb);
	fputs(prea, fp_pa);
	fputs(bait, fp_bb);
	char *sm = parse_rg(hdr, "SM");
	char *lb = parse_rg(hdr, "LB");
	char c[4] = {0}, r, a;
	for (int i = 0; i < 4; ++i)
	{
		r = base[i];
		for (int j = 0; j < 4; ++j)
		{
			if (i == j) continue;
			a = base[j];
			for (int k = 0; k < 4; ++k)
			{
				for (int l = 0; l < 4; ++l)
				{
					c[0] = base[k];
					c[1] = base[i];
					c[2] = base[l];
					int pr = m[1][seq2idx(r, c)] + m[2][seq2idx(r, c)] + m[0][qes2idx(r, c)] + m[3][qes2idx(r, c)];
					int pa = m[1][seq2idx(a, c)] + m[2][seq2idx(a, c)] + m[0][qes2idx(a, c)] + m[3][qes2idx(a, c)];
					int cr = m[0][seq2idx(r, c)] + m[3][seq2idx(r, c)] + m[1][qes2idx(r, c)] + m[2][qes2idx(r, c)];
					int ca = m[0][seq2idx(a, c)] + m[3][seq2idx(a, c)] + m[1][qes2idx(a, c)] + m[2][qes2idx(a, c)];
					int sum = pr + pa + cr + ca;
					double er = (sum == 0) ? pow(10, -10) : max(pow(10, -10), (double)(pa - ca)/sum);
					int qs = (int)round(-10*log10(er));
					fprintf(fp_pa, "%s\t%s\t%c\t%c\t%s\t%d\t%d\t%d\t%d\t%.*f\t%d\n", sm, lb, r, a, c, pr, pa, cr, ca, er > pow(10, -10) ? 6 : 0, er, qs);
					int fr = m[0][seq2idx(r, c)] + m[1][seq2idx(r, c)] + m[2][seq2idx(r, c)] + m[3][seq2idx(r, c)];
					int fa = m[0][seq2idx(a, c)] + m[1][seq2idx(a, c)] + m[2][seq2idx(a, c)] + m[3][seq2idx(a, c)];
					int rr = m[0][qes2idx(r, c)] + m[1][qes2idx(r, c)] + m[2][qes2idx(r, c)] + m[3][qes2idx(r, c)];
					int ra = m[0][qes2idx(a, c)] + m[1][qes2idx(a, c)] + m[2][qes2idx(a, c)] + m[3][qes2idx(a, c)];
					double f_er = (fr + fa) > 0 ? max(pow(10, -10), (double)fa/(fa+fr)) : pow(10, -10);
					double r_er = (rr + ra) > 0 ? max(pow(10, -10), (double)ra/(ra+rr)) : pow(10, -10);
					er = max(pow(10, -10), f_er - r_er);
					qs = (int)round(-10*log10(er));
					fprintf(fp_bb, "%s\t%s\t%c\t%c\t%s\t%d\t%d\t%d\t%d\t%.*f\t%.*f\t%.*f\t%d\n", sm, lb, r, a, c, fr, fa, rr, ra, f_er > pow(10, -10) ? 6 : 0, f_er, r_er > pow(10, -10) ? 6 : 0, r_er, er > pow(10, -10) ? 6 : 0, er, qs);
				}
			}
		}
	}
	free(sm);
	free(lb);
	free(ks_pa.s);
	free(ks_bb.s);
	fclose(fp_pa);
	fclose(fp_bb);
}

char *get_read(const bam1_t *b)
{
	int len = b->core.l_qseq + 1;
	char *read = calloc(1, len);
	char *seq = (char *)bam_get_seq(b);
	if (!read) return NULL;
	for (int n = 0; n < b->core.l_qseq; ++n)
		read[n] = seq_nt16_str[bam_seqi(seq,n)];
	return read;
}

int get_qual(const bam1_t *b, char **qual_out)
{
	char *quality = calloc(1, b->core.l_qseq + 1);
	char *q = (char *)bam_get_qual(b);
	if (!quality) return -1;
	if (*q == '\xff')
	{
		free(quality);
		*qual_out = NULL;
		return 0;
	}
	for (int n = 0; n < b->core.l_qseq; ++n)
		quality[n] = q[n] + 33;
	*qual_out = quality;
	return 0;
}

int seq2idx(const char q, const char *r)
{
	// seq_nt16_table  ATGC -> 1248
	// seq_nt16_int    1248 -> 0123
	uint8_t idx = 0;
	idx |= (uint8_t)seq_nt16_int[seq_nt16_table[q]] << 6;
	for (int i = 0; i < 3; ++i)
		idx |= (uint8_t)(seq_nt16_int[seq_nt16_table[r[i]]] << (2-i)*2);
	return idx;
}

int qes2idx(const char q, const char *r)
{
	uint8_t idx = 0;
	idx |= (uint8_t)(~seq_nt16_int[seq_nt16_table[q]] & 0x3) << 6;
	for (int i = 0; i < 3; ++i)
		idx |= (uint8_t)((~seq_nt16_int[seq_nt16_table[r[i]]] & 0x3) << i*2);
	return idx;
}

void idx2seq(int n, char *q, char *r)
{
	q[0] = base[(uint8_t)((n >> 6) & 0x3)];
	for (int i = 0; i < 3; ++i)
		r[i] = base[(uint8_t)((n >> (2-i)*2) & 0x3)];
}

char *get_cigar(bam1_t *b)
{
	kstring_t cigar = {0, 0, 0};
	uint32_t *ez = bam1_cigar(b);
	int32_t on = b->core.n_cigar;
	for (int i = 0; i < on; ++i)
	{
		int ol = bam_cigar_oplen(ez[i]), op = bam_cigar_op(ez[i]);
		ksprintf(&cigar, "%d%c", ol, "MID"[(uint8_t)op]);
	}
	return(ks_release(&cigar));
}

static inline void nt_to_upper(char *nt)
{
    if ( *nt >= 97 ) *nt -= 32;
}

uint64_t tot_reads(const char *fn)
{
	uint64_t tt = 0;
	BGZF *in = bgzf_open(fn, "r"); assert(in);
	bam_hdr_t *hdr = bam_hdr_read(in);
	bam_index_t *idx = bam_index_load(fn);
	for (int i = 0; i < hdr->n_targets; ++i)
	{
		uint64_t u, v;
		hts_idx_get_stat(idx, i, &u, &v);
		tt += u + v;
	}
	tt += hts_idx_get_n_no_coor(idx);
	bam_hdr_destroy(hdr);
	bam_index_destroy(idx);
	bgzf_close(in);
	return tt;
}

void progressbar(uint64_t _c, uint64_t _t, int _w)
{
	int c = _w / 3 * _c / _t;
	int t = _w / 3;
	if (isatty(fileno(stdin)))
	{
		wchar_t bar[(t + 1) * sizeof(wchar_t)];
		for (int i = 0; i < c; ++i)
			bar[i] = (wchar_t)0x25FC;
		for (int i = c; i < t; ++i)
			bar[i] = (wchar_t)0x25FB;
		bar[t] = '\0';
		fprintf(stderr, "\r\e[2m%ls%3d%%\e%s[0m", bar, min(c*100/t, 100), c == t ? "\n" : "");
		fflush(stderr);
	}
}

char *parse_rg(const bam_hdr_t *hdr, const char *tag)
{
	char *s;
	kstring_t ks = {0, 0, 0}, tg = {0, 0, 0};
	ksprintf(&tg, "\t%s:", tag);
	char *smp = strstr(hdr->text, tg.s);
	free(tg.s);
	if (!smp)
		s = strdup("N/A");
	else
	{
		smp += 4;
		while (*smp != '\t' && *smp != ' ' && *smp != '\n')
			kputc(*smp++, &ks);
		s = ks_release(&ks);
	}
	return s;
}

void ini_arg(arg_t *arg)
{
	arg->bam = 0;
	arg->ref = 0;
	arg->pfx = 0;
	arg->max_is = 600;
	arg->min_is = 60;
	arg->min_mq = 30;
	arg->min_qs = 20;
}

void usage(char *str)
{
	fprintf(stderr, "\nCollect Sequencing Artifact Metrics\n");
	fprintf(stderr, "\nUsage:  \e[1;31m%s\e[0;0m -i <bam> -r <ref> -p <prefix>\n", basename(str));
	fprintf(stderr, "\n Required:\n");
	fprintf(stderr, "    -i --bam    [STR] Input bam file (bai index required)\n");
	fprintf(stderr, "    -r --ref    [STR] reference fa (fai index required)\n");
	fprintf(stderr, "\n Optional:\n");
	fprintf(stderr, "    -p --pfx    [STR] output prefix (default \"sample\")\n");
	fprintf(stderr, "    -M --max_is [INT] max insert size (default 600)\n");
	fprintf(stderr, "    -m --min_is [INT] min insert size (default 60)\n");
	fprintf(stderr, "    -q --min_mq [INT] min mapping quality (default 30)\n");
	fprintf(stderr, "    -Q --min_qs [INT] min base quality score (default 20)\n");
	fputc('\n', stderr);
	fprintf(stderr, "    -v --version      show version\n");
	fprintf(stderr, "    -h --help         show help message\n");
	fprintf(stderr, "\nContact:  slw287r@163.com\n\n");
	exit(EXIT_FAILURE);
}
