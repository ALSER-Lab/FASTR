/* The MIT License

   Copyright (c) 2008, 2009, 2011 Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Last Modified: 05MAR2012 */

#ifndef AC_KSEQ_H
#define AC_KSEQ_H

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
///CHANGE: ALSER LAB///
#include <stdint.h>
///END OF CHANGE///

#define KS_SEP_SPACE 0 // isspace(): \t, \n, \v, \f, \r
#define KS_SEP_TAB   1 // isspace() && !' '
#define KS_SEP_LINE  2 // line separator: "\n" (Unix) or "\r\n" (Windows)
#define KS_SEP_MAX   2

#define __KS_TYPE(type_t)						\
	typedef struct __kstream_t {				\
		unsigned char *buf;						\
		int begin, end, is_eof;					\
		type_t f;								\
	} kstream_t;

#define ks_eof(ks) ((ks)->is_eof && (ks)->begin >= (ks)->end)
#define ks_rewind(ks) ((ks)->is_eof = (ks)->begin = (ks)->end = 0)

#define __KS_BASIC(type_t, __bufsize)								\
	static inline kstream_t *ks_init(type_t f)						\
	{																\
		kstream_t *ks = (kstream_t*)calloc(1, sizeof(kstream_t));	\
		ks->f = f;													\
		ks->buf = (unsigned char*)malloc(__bufsize);				\
		return ks;													\
	}																\
	static inline void ks_destroy(kstream_t *ks)					\
	{																\
		if (ks) {													\
			free(ks->buf);											\
			free(ks);												\
		}															\
	}

#define __KS_GETC(__read, __bufsize)						\
	static inline int ks_getc(kstream_t *ks)				\
	{														\
		if (ks->is_eof && ks->begin >= ks->end) return -1;	\
		if (ks->begin >= ks->end) {							\
			ks->begin = 0;									\
			ks->end = __read(ks->f, ks->buf, __bufsize);	\
			if (ks->end == 0) { ks->is_eof = 1; return -1;}	\
		}													\
		return (int)ks->buf[ks->begin++];					\
	}

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
  size_t l, m;
  char *s;
} kstring_t;
#endif

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#define __KS_GETUNTIL(__read, __bufsize)								\
	static int ks_getuntil2(kstream_t *ks, int delimiter, kstring_t *str, int *dret, int append) \
	{																	\
		int gotany = 0;													\
		if (dret) *dret = 0;											\
		str->l = append? str->l : 0;									\
		for (;;) {														\
			int i;														\
			if (ks->begin >= ks->end) {									\
				if (!ks->is_eof) {										\
					ks->begin = 0;										\
					ks->end = __read(ks->f, ks->buf, __bufsize);		\
					if (ks->end == 0) { ks->is_eof = 1; break; }		\
				} else break;											\
			}															\
			if (delimiter == KS_SEP_LINE) { \
				for (i = ks->begin; i < ks->end; ++i) \
					if (ks->buf[i] == '\n') break; \
			} else if (delimiter > KS_SEP_MAX) {						\
				for (i = ks->begin; i < ks->end; ++i)					\
					if (ks->buf[i] == delimiter) break;					\
			} else if (delimiter == KS_SEP_SPACE) {						\
				for (i = ks->begin; i < ks->end; ++i)					\
					if (isspace(ks->buf[i])) break;						\
			} else if (delimiter == KS_SEP_TAB) {						\
				for (i = ks->begin; i < ks->end; ++i)					\
					if (isspace(ks->buf[i]) && ks->buf[i] != ' ') break; \
			} else i = 0; /* never come to here! */						\
			if (str->m - str->l < (size_t)(i - ks->begin + 1)) {		\
				str->m = str->l + (i - ks->begin) + 1;					\
				kroundup32(str->m);										\
				str->s = (char*)realloc(str->s, str->m);				\
			}															\
			gotany = 1;													\
			memcpy(str->s + str->l, ks->buf + ks->begin, i - ks->begin); \
			str->l = str->l + (i - ks->begin);							\
			ks->begin = i + 1;											\
			if (i < ks->end) {											\
				if (dret) *dret = ks->buf[i];							\
				break;													\
			}															\
		}																\
		if (!gotany && ks_eof(ks)) return -1;							\
		if (str->s == 0) {												\
			str->m = 1;													\
			str->s = (char*)calloc(1, 1);								\
		} else if (delimiter == KS_SEP_LINE && str->l > 1 && str->s[str->l-1] == '\r') --str->l; \
		str->s[str->l] = '\0';											\
		return str->l;													\
	} \
	static inline int ks_getuntil(kstream_t *ks, int delimiter, kstring_t *str, int *dret) \
	{ return ks_getuntil2(ks, delimiter, str, dret, 0); }

#define KSTREAM_INIT(type_t, __read, __bufsize) \
	__KS_TYPE(type_t)							\
	__KS_BASIC(type_t, __bufsize)				\
	__KS_GETC(__read, __bufsize)				\
	__KS_GETUNTIL(__read, __bufsize)

#define kseq_rewind(ks) ((ks)->last_char = (ks)->f->is_eof = (ks)->f->begin = (ks)->f->end = 0)

///CHANGE: ALSER LAB///
#ifndef klib_unused
#if (defined __clang__ && __clang_major__ >= 3) || (defined __GNUC__ && __GNUC__ >= 3)
#define klib_unused __attribute__ ((__unused__))
#else
#define klib_unused
#endif
#endif

static inline unsigned char kseq_fastr_unescape(unsigned char b)
{
	return b == 255u ? 10u : b == 254u ? 64u : b;
}

static klib_unused int kseq_fastr__parse_int_csv(const char *s, int *out, int maxn)
{
	int n = 0; const char *p = s;
	while (*p && n < maxn) {
		while (*p == ' ' || *p == ',' || *p == '\t') ++p;
		if (!*p) break;
		char *e = NULL; long v = strtol(p, &e, 10);
		if (e == p) break;
		out[n++] = (int)v; p = e;
	}
	return n;
}

static klib_unused void kseq_fastr__build_inverse(int32_t *inverse64, const char *qmap_csv)
{
	uint8_t lut[94]; size_t n = 0; const char *p = qmap_csv;
	uint8_t max_slot = 0; size_t k;
	int i; for (i = 0; i < 64; ++i) inverse64[i] = 0;
	if (!qmap_csv || !*qmap_csv) return;
	while (*p && n < 94) {
		while (*p == ' ' || *p == ',' || *p == '\t') ++p;
		if (!*p) break;
		char *e = NULL; long v = strtol(p, &e, 10);
		if (e == p) break;
		if (v < 0) v = 0;
		if (v > 63) v = 63;
		lut[n++] = (uint8_t)v; p = e;
	}
	for (k = 0; k < n; ++k) if (lut[k] > max_slot) max_slot = lut[k];
	for (k = 0; k < n; ++k) { uint8_t s = lut[k]; if (s >= 1 && s <= 63) inverse64[s-1] = (int32_t)k; }
	if (max_slot >= 1 && max_slot <= 63) {
		int32_t si = max_slot - 1; int32_t rep = -1;
		for (k = 0; k < n; ++k) if (lut[k] == max_slot) { rep = (int32_t)k; break; }
		inverse64[si] = rep < 0 ? 0 : rep;
	}
}

static klib_unused void kseq_fastr__apply_decode(int32_t *inverse64, const char *qdec_csv)
{
	size_t i = 0; const char *p = qdec_csv;
	if (!qdec_csv || !*qdec_csv) return;
	while (*p && i < 64) {
		while (*p == ' ' || *p == '\t' || *p == ',') ++p;
		if (!*p) break;
		char *e = NULL; long v = strtol(p, &e, 10);
		if (e == p) break;
		if (v < 0) v = 0;
		if (v > 93) v = 93;
		inverse64[i++] = (int32_t)v; p = e;
	}
}

static klib_unused void kseq_fastr__init_nonet_addtab(int addtab[128][7])
{
	int m, j;
	for (m = 0; m < 128; ++m) for (j = 0; j < 7; ++j) addtab[m][j] = 63 * ((m >> j) & 1);
}
///END OF CHANGE///

#define __KSEQ_BASIC(SCOPE, type_t)										\
	SCOPE kseq_t *kseq_init(type_t fd)									\
	{																	\
		kseq_t *s = (kseq_t*)calloc(1, sizeof(kseq_t));					\
		s->f = ks_init(fd);												\
		return s;														\
	}																	\
	SCOPE void kseq_destroy(kseq_t *ks)									\
	{																	\
		if (!ks) return;												\
		free(ks->name.s); free(ks->comment.s); free(ks->seq.s);	free(ks->qual.s); \
		/* ///CHANGE: ALSER LAB/// */ \
		free(ks->fastr_raw_scratch.s); \
		/* ///END OF CHANGE/// */ \
		ks_destroy(ks->f);												\
		free(ks);														\
	}

///CHANGE: ALSER LAB///
#define __KSEQ_FASTR_ENABLE(SCOPE) \
	SCOPE void kseq_fastr__parse_header(kseq_t *ks) \
	{ \
		int gray[5] = {0, 1, 64, 127, 190}; \
		kstring_t line = {0, 0, 0}; \
		char *qmap = NULL; \
		char *qdec = NULL; \
		int i, dret; \
		ks->fastr_enc = 0; \
		ks->fastr_max_phred = 93; \
		ks->fastr_n_phred = 2; \
		for (i = 0; ; ++i) { \
			if (i > 0) { \
				int nc = ks_getc(ks->f); \
				if (nc != '#') { if (nc >= 0) --ks->f->begin; break; } \
			} \
			if (ks_getuntil(ks->f, KS_SEP_LINE, &line, &dret) < 0) break; \
			if (!line.s) continue; \
			if (strncmp(line.s, "ENCODING=", 9) == 0 || strncmp(line.s, "#ENCODING=", 10) == 0) { \
				const char *v = strchr(line.s, '=') + 1; \
				if      (strncmp(v, "nibble", 6) == 0) ks->fastr_enc = 1; \
				else if (strncmp(v, "nonet", 5) == 0)  ks->fastr_enc = 2; \
				else                                    ks->fastr_enc = 0; \
			} else if (strncmp(line.s, "GRAY_VALS=", 10) == 0 || strncmp(line.s, "#GRAY_VALS=", 11) == 0) { \
				kseq_fastr__parse_int_csv(strchr(line.s, '=') + 1, gray, 5); \
			} else if (strncmp(line.s, "QUALITY_MAP=", 12) == 0 || strncmp(line.s, "#QUALITY_MAP=", 13) == 0) { \
				free(qmap); qmap = strdup(strchr(line.s, '=') + 1); \
			} else if (strncmp(line.s, "QUALITY_DECODE=", 15) == 0 || strncmp(line.s, "#QUALITY_DECODE=", 16) == 0) { \
				free(qdec); qdec = strdup(strchr(line.s, '=') + 1); \
			} else if (strncmp(line.s, "N_QUALITY=", 10) == 0 || strncmp(line.s, "#N_QUALITY=", 11) == 0) { \
				char *v = strchr(line.s, '='); \
				if (v && v[1]) ks->fastr_n_phred = (int)v[1] - 33; \
			} \
		} \
		if (line.s) free(line.s); \
		{ \
			int gA=gray[1], gC=gray[2], gG=gray[3], gT=gray[4], j; \
			unsigned char base[256], band_start[256]; \
			for (j = 0; j < 256; ++j) { base[j] = 'N'; band_start[j] = 0; } \
			for (j = gA; j < gC && j < 256; ++j) { base[j] = 'A'; band_start[j] = (unsigned char)gA; } \
			for (j = gC; j < gG && j < 256; ++j) { base[j] = 'C'; band_start[j] = (unsigned char)gC; } \
			for (j = gG; j < gT && j < 256; ++j) { base[j] = 'G'; band_start[j] = (unsigned char)gG; } \
			for (j = gT; j < gT + 63 && j < 253; ++j) { base[j] = 'T'; band_start[j] = (unsigned char)gT; } \
			for (j = 0; j < 256; ++j) { \
				unsigned char ub = kseq_fastr_unescape((unsigned char)j); \
				base[j] = base[ub]; \
				band_start[j] = band_start[ub]; \
			} \
			kseq_fastr__build_inverse(ks->fastr_inverse, qmap); \
			kseq_fastr__apply_decode(ks->fastr_inverse, qdec); \
			free(qmap); \
			free(qdec); \
			for (j = 0; j < 256; ++j) { \
				int v; unsigned char c = base[j]; \
				if (c == 'N') v = ks->fastr_n_phred; \
				else { \
					int y = (int)kseq_fastr_unescape((unsigned char)j) - (int)band_start[j]; \
					if (y < 0) y = 0; \
					if (y > 63) y = 63; \
					v = ks->fastr_inverse[y > 62 ? 62 : y]; \
				} \
				if (v < 0) v = 0; \
				if (v > ks->fastr_max_phred) v = ks->fastr_max_phred; \
				ks->fastr_byte_lut[j] = (uint16_t)((unsigned char)c | ((unsigned char)(33 + v) << 8)); \
				{ \
					int q; \
					if (c == 'N') q = ks->fastr_n_phred; \
					else { q = (int)kseq_fastr_unescape((unsigned char)j) - (int)band_start[j]; if (q < 0) q = 0; if (q > ks->fastr_max_phred) q = ks->fastr_max_phred; } \
					ks->fastr_nonet_lut[j] = (uint16_t)((unsigned char)c | ((unsigned char)q << 8)); \
				} \
			} \
			{ \
				unsigned char nib_base[16]; int nib_slot[16]; \
				for (j = 0; j < 16; ++j) { nib_base[j] = 0; nib_slot[j] = 0; } \
				nib_base[0] = 'N'; \
				for (j = 1; j <= 3; ++j)  { nib_base[j] = 'A'; nib_slot[j] = j - 1; } \
				for (j = 4; j <= 6; ++j)  { nib_base[j] = 'C'; nib_slot[j] = j - 4; } \
				for (j = 7; j <= 9; ++j)  { nib_base[j] = 'G'; nib_slot[j] = j - 7; } \
				for (j = 10; j <= 12; ++j) { nib_base[j] = 'T'; nib_slot[j] = j - 10; } \
				for (j = 0; j < 256; ++j) { \
					unsigned char ub = kseq_fastr_unescape((unsigned char)j); \
					unsigned char hi = (unsigned char)(ub >> 4), lo = (unsigned char)(ub & 0x0F); \
					uint32_t packed = 0; \
					if (hi <= 12 && nib_base[hi]) { \
						int v = (nib_base[hi]=='N') ? ks->fastr_n_phred : ks->fastr_inverse[nib_slot[hi]]; \
						if (v < 0) v = 0; \
						if (v > ks->fastr_max_phred) v = ks->fastr_max_phred; \
						packed |= (uint32_t)(unsigned char)nib_base[hi] | ((uint32_t)(unsigned char)(33+v) << 8); \
					} \
					if (lo <= 12 && nib_base[lo]) { \
						int v = (nib_base[lo]=='N') ? ks->fastr_n_phred : ks->fastr_inverse[nib_slot[lo]]; \
						if (v < 0) v = 0; \
						if (v > ks->fastr_max_phred) v = ks->fastr_max_phred; \
						packed |= (uint32_t)(unsigned char)nib_base[lo] << 16 | ((uint32_t)(unsigned char)(33+v) << 24); \
					} \
					ks->fastr_nib_lut[j] = packed; \
					ks->fastr_nib_base2[j] = (uint16_t)((packed & 0xFFu) | (((packed >> 16) & 0xFFu) << 8)); \
					ks->fastr_nib_qual2[j] = (uint16_t)(((packed >> 8) & 0xFFu) | (((packed >> 24) & 0xFFu) << 8)); \
					ks->fastr_nib_cnt[j] = (unsigned char)(((packed & 0xFFu) != 0) + ((packed & 0xFF0000u) != 0)); \
				} \
			} \
			kseq_fastr__init_nonet_addtab(ks->fastr_nonet_addtab); \
		} \
	}

#define __KSEQ_FASTR_READ(SCOPE) \
	SCOPE int kseq_fastr__read_record(kseq_t *seq) \
	{ \
		kstream_t *ks = seq->f; \
		int dret; \
		seq->comment.l = seq->seq.l = seq->qual.l = 0; \
		if (ks_getc(ks) < 0) return -1; \
		if (ks_getuntil2(ks, '\n', &seq->name, &dret, 0) < 0) return -1; \
		{ char *sp=(char*)memchr(seq->name.s,' ',seq->name.l); \
			if (sp){ size_t nl=(size_t)(sp-seq->name.s),cl=seq->name.l-nl-1; \
				if(seq->comment.m<cl+1){seq->comment.m=cl+1;kroundup32(seq->comment.m);seq->comment.s=(char*)realloc(seq->comment.s,seq->comment.m);} \
				memcpy(seq->comment.s,sp+1,cl);seq->comment.s[cl]='\0';seq->comment.l=cl; \
				seq->name.s[nl]='\0';seq->name.l=nl; } } \
		{ \
			kstring_t *raw = &seq->fastr_raw_scratch; \
			size_t L, i; \
			if (ks_getuntil2(ks, '\n', raw, 0, 0) < 0) return -2; \
			{ \
				size_t need = (seq->fastr_enc == 1 ? raw->l * 2 : raw->l) + 2; \
				if (seq->seq.m < need) { seq->seq.m = need; kroundup32(seq->seq.m); seq->seq.s = (char*)realloc(seq->seq.s, seq->seq.m); } \
				if (seq->qual.m < need) { seq->qual.m = need; kroundup32(seq->qual.m); seq->qual.s = (char*)realloc(seq->qual.s, seq->qual.m); } \
			} \
			if (seq->fastr_enc == 1) { \
				size_t nfull = raw->l ? raw->l - 1 : 0; \
				uint32_t v; \
				L = 0; \
				for (i = 0; i < nfull; ++i) { \
					unsigned char b = (unsigned char)raw->s[i]; \
					*(uint16_t*)(seq->seq.s + L) = seq->fastr_nib_base2[b]; \
					*(uint16_t*)(seq->qual.s + L) = seq->fastr_nib_qual2[b]; \
					L += seq->fastr_nib_cnt[b]; \
				} \
				if (raw->l) { \
					v = seq->fastr_nib_lut[(unsigned char)raw->s[raw->l - 1]]; \
					if (v & 0xFFu)     { seq->seq.s[L] = (char)(v & 0xFF);         seq->qual.s[L] = (char)((v >> 8) & 0xFF);  ++L; } \
					if (v & 0xFF0000u) { seq->seq.s[L] = (char)((v >> 16) & 0xFF); seq->qual.s[L] = (char)((v >> 24) & 0xFF); ++L; } \
				} \
			} else if (seq->fastr_enc == 2) { \
				unsigned char *p = (unsigned char*)raw->s; \
				unsigned char *sep = (unsigned char*)memchr(p, 0xFD, raw->l); \
				L = sep ? (size_t)(sep - p) : raw->l; \
				const unsigned char *bit = sep ? sep + 1 : NULL; \
				const int mp = seq->fastr_max_phred; \
				if (bit) { \
					size_t k = 0; \
					while (k + 7 <= L) { \
						const int *add = seq->fastr_nonet_addtab[(*bit++) & 0x7F]; \
						int j; for (j = 0; j < 7; ++j) { \
							uint16_t v = seq->fastr_nonet_lut[p[k+j]]; \
							int ph = (int)(v >> 8) + add[j]; if (ph > mp) ph = mp; \
							seq->seq.s[k+j] = (char)(v & 0xFF); \
							seq->qual.s[k+j] = (char)(33 + ph); \
						} \
						k += 7; \
					} \
					if (k < L) { unsigned mmask = *bit; int j = 0; \
						for (; k < L; ++k, ++j) { \
							uint16_t v = seq->fastr_nonet_lut[p[k]]; \
							int ph = (int)(v >> 8) + 63 * ((mmask >> j) & 1u); if (ph > mp) ph = mp; \
							seq->seq.s[k] = (char)(v & 0xFF); \
							seq->qual.s[k] = (char)(33 + ph); \
						} \
					} \
				} else { \
					for (i = 0; i < L; ++i) { \
						uint16_t v = seq->fastr_nonet_lut[p[i]]; \
						seq->seq.s[i] = (char)(v & 0xFF); \
						seq->qual.s[i] = (char)(33 + (v >> 8)); \
					} \
				} \
			} else { \
				L = raw->l; \
				for (i = 0; i < L; ++i) { \
					uint16_t v = seq->fastr_byte_lut[(unsigned char)raw->s[i]]; \
					seq->seq.s[i] = (char)(v & 0xFF); \
					seq->qual.s[i] = (char)(v >> 8); \
				} \
			} \
			seq->seq.l = seq->qual.l = L; \
			seq->seq.s[L] = '\0'; seq->qual.s[L] = '\0'; \
		} \
		return (int)seq->seq.l; \
	}
///END OF CHANGE///

/* Return value:
   >=0  length of the sequence (normal)
   -1   end-of-file
   -2   truncated quality string
 */
#define __KSEQ_READ(SCOPE) \
	SCOPE int kseq_read(kseq_t *seq) \
	{ \
		int c; \
		kstream_t *ks = seq->f; \
		/* ///CHANGE: ALSER LAB/// */ \
		if (!seq->fastr_checked) { \
			seq->fastr_checked = 1; \
			c = ks_getc(ks); \
			if (c == '#') { \
				seq->is_fastr = 1; \
				kseq_fastr__parse_header(seq); \
				return kseq_fastr__read_record(seq); \
			} else if (c == -1) { \
				return -1; \
			} else if (c == '>' || c == '@') { \
				seq->last_char = c; \
			} else { \
				while ((c = ks_getc(ks)) != -1 && c != '>' && c != '@'); \
				if (c == -1) return -1; \
				seq->last_char = c; \
			} \
		} else if (seq->is_fastr) { \
			return kseq_fastr__read_record(seq); \
		} \
		/* ///END OF CHANGE/// */ \
		if (seq->last_char == 0) { /* then jump to the next header line */ \
			while ((c = ks_getc(ks)) != -1 && c != '>' && c != '@'); \
			if (c == -1) return -1; /* end of file */ \
			seq->last_char = c; \
		} /* else: the first header char has been read in the previous call */ \
		seq->comment.l = seq->seq.l = seq->qual.l = 0; /* reset all members */ \
		if (ks_getuntil(ks, 0, &seq->name, &c) < 0) return -1; /* normal exit: EOF */ \
		if (c != '\n') ks_getuntil(ks, KS_SEP_LINE, &seq->comment, 0); /* read FASTA/Q comment */ \
		if (seq->seq.s == 0) { /* we can do this in the loop below, but that is slower */ \
			seq->seq.m = 256; \
			seq->seq.s = (char*)malloc(seq->seq.m); \
		} \
		while ((c = ks_getc(ks)) != -1 && c != '>' && c != '+' && c != '@') { \
			if (c == '\n') continue; /* skip empty lines */ \
			seq->seq.s[seq->seq.l++] = c; /* this is safe: we always have enough space for 1 char */ \
			ks_getuntil2(ks, KS_SEP_LINE, &seq->seq, 0, 1); /* read the rest of the line */ \
		} \
		if (c == '>' || c == '@') seq->last_char = c; /* the first header char has been read */	\
		if (seq->seq.l + 1 >= seq->seq.m) { /* seq->seq.s[seq->seq.l] below may be out of boundary */ \
			seq->seq.m = seq->seq.l + 2; \
			kroundup32(seq->seq.m); /* rounded to the next closest 2^k */ \
			seq->seq.s = (char*)realloc(seq->seq.s, seq->seq.m); \
		} \
		seq->seq.s[seq->seq.l] = 0;	/* null terminated string */ \
		if (c != '+') return seq->seq.l; /* FASTA */ \
		if (seq->qual.m < seq->seq.m) {	/* allocate memory for qual in case insufficient */ \
			seq->qual.m = seq->seq.m; \
			seq->qual.s = (char*)realloc(seq->qual.s, seq->qual.m); \
		} \
		while ((c = ks_getc(ks)) != -1 && c != '\n'); /* skip the rest of '+' line */ \
		if (c == -1) return -2; /* error: no quality string */ \
		while (ks_getuntil2(ks, KS_SEP_LINE, &seq->qual, 0, 1) >= 0 && seq->qual.l < seq->seq.l); \
		seq->last_char = 0;	/* we have not come to the next header line */ \
		if (seq->seq.l != seq->qual.l) return -2; /* error: qual string is of a different length */ \
		return seq->seq.l; \
	}

#define __KSEQ_TYPE(type_t)						\
	typedef struct {							\
		kstring_t name, comment, seq, qual;		\
		int last_char;							\
		kstream_t *f;							\
		/* ///CHANGE: ALSER LAB/// */ \
		int fastr_checked, is_fastr, fastr_enc; \
		int fastr_max_phred, fastr_n_phred; \
		kstring_t fastr_raw_scratch; \
		int32_t fastr_inverse[64]; \
		uint16_t fastr_byte_lut[256]; \
		uint16_t fastr_nonet_lut[256]; \
		uint32_t fastr_nib_lut[256]; \
		uint16_t fastr_nib_base2[256]; \
		uint16_t fastr_nib_qual2[256]; \
		unsigned char fastr_nib_cnt[256]; \
		int fastr_nonet_addtab[128][7]; \
		/* ///END OF CHANGE/// */ \
	} kseq_t;

#define KSEQ_INIT2(SCOPE, type_t, __read)		\
	KSTREAM_INIT(type_t, __read, 16384)			\
	__KSEQ_TYPE(type_t)							\
	__KSEQ_BASIC(SCOPE, type_t)					\
	/* ///CHANGE: ALSER LAB/// */ \
	__KSEQ_FASTR_ENABLE(SCOPE) \
	__KSEQ_FASTR_READ(SCOPE) \
	/* ///END OF CHANGE/// */ \
	__KSEQ_READ(SCOPE)

#define KSEQ_INIT(type_t, __read) KSEQ_INIT2(static, type_t, __read)

#define KSEQ_DECLARE(type_t) \
	__KS_TYPE(type_t) \
	__KSEQ_TYPE(type_t) \
	extern kseq_t *kseq_init(type_t fd); \
	void kseq_destroy(kseq_t *ks); \
	int kseq_read(kseq_t *seq);

#endif
