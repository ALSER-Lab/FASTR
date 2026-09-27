///CHANGE: ALSER LAB///
package uk.ac.babraham.FastQC.Sequence;

import java.io.*;
import java.nio.charset.StandardCharsets;

import uk.ac.babraham.FastQC.FastQCConfig;

public class FASTRFile implements SequenceFile {

    private static final int BUFSIZE = 1 << 16;

    private final File file;
    private final String name;
    private final FileInputStream fis;
    private final long fileSize;

    private final byte[] buf = new byte[BUFSIZE];
    private int begin = 0, end = 0;
    private boolean isEof = false;

    private byte[] line = new byte[1 << 12];
    private int lineLen = 0;
    private byte[] seqBuf = new byte[1 << 12];
    private byte[] qualBuf = new byte[1 << 12];

    private int enc = 0;
    private int maxPhred = 93;
    private int nPhred = 2;
    private final int[] inverse = new int[64];
    private final short[] byteLut = new short[256];
    private final short[] nonetLut = new short[256];
    private final short[] nibBase2 = new short[256];
    private final short[] nibQual2 = new short[256];
    private final byte[] nibCnt = new byte[256];
    private final int[] nibLut = new int[256];

    private Sequence nextSequence;

    protected FASTRFile(FastQCConfig config, File file) throws SequenceFormatException, IOException {
        this.file = file;
        this.name = file.getName();
        this.fileSize = file.length();
        this.fis = new FileInputStream(file);

        int c = getc();
        if (c != '#') throw new SequenceFormatException("Not a FASTR file (first byte is not '#'): " + name);
        parseHeader();
        readNext();
    }

    private int fill() throws IOException {
        begin = 0;
        int n = 0;
        while (n < BUFSIZE) {
            int r = fis.read(buf, n, BUFSIZE - n);
            if (r < 0) { isEof = true; break; }
            n += r;
        }
        end = n;
        return n;
    }

    private int getc() throws IOException {
        if (begin >= end) {
            if (isEof) return -1;
            if (fill() == 0) return -1;
        }
        return buf[begin++] & 0xFF;
    }

    private int getLine() throws IOException {
        lineLen = 0;
        if (begin >= end && isEof) return -1;
        boolean gotAny = false;
        for (;;) {
            if (begin >= end) {
                if (isEof || fill() == 0) break;
            }
            int i = begin;
            while (i < end && buf[i] != '\n') ++i;
            int n = i - begin;
            if (lineLen + n > line.length) {
                byte[] nl = new byte[Math.max(line.length * 2, lineLen + n)];
                System.arraycopy(line, 0, nl, 0, lineLen);
                line = nl;
            }
            System.arraycopy(buf, begin, line, lineLen, n);
            lineLen += n;
            gotAny = true;
            begin = i + 1;
            if (i < end) return lineLen;
        }
        return gotAny ? lineLen : -1;
    }

    private static int unescape(int b) {
        return b == 255 ? 10 : b == 254 ? 64 : b;
    }

    private static int[] parseIntCsv(String s, int maxn) {
        int[] out = new int[maxn];
        int n = 0;
        for (String t : s.split("[,\\s]+")) {
            if (t.isEmpty()) continue;
            if (n >= maxn) break;
            try { out[n++] = Integer.parseInt(t); } catch (NumberFormatException e) { break; }
        }
        int[] r = new int[n];
        System.arraycopy(out, 0, r, 0, n);
        return r;
    }

    private static String value(String l, String key) {
        String k = key + "=";
        if (l.startsWith(k)) return l.substring(k.length());
        if (l.startsWith("#" + k)) return l.substring(k.length() + 1);
        return null;
    }

    private void parseHeader() throws IOException {
        int[] gray = {0, 1, 64, 127, 190};
        String qmap = null, qdec = null, v;
        for (int i = 0; ; ++i) {
            if (i > 0) {
                int nc = getc();
                if (nc != '#') { if (nc >= 0) --begin; break; }
            }
            if (getLine() < 0) break;
            String l = new String(line, 0, lineLen, StandardCharsets.ISO_8859_1);
            if ((v = value(l, "ENCODING")) != null) {
                enc = v.startsWith("nibble") ? 1 : v.startsWith("nonet") ? 2 : 0;
            } else if ((v = value(l, "GRAY_VALS")) != null) {
                int[] g = parseIntCsv(v, 5);
                System.arraycopy(g, 0, gray, 0, g.length);
            } else if ((v = value(l, "QUALITY_MAP")) != null) {
                qmap = v;
            } else if ((v = value(l, "QUALITY_DECODE")) != null) {
                qdec = v;
            } else if ((v = value(l, "N_QUALITY")) != null) {
                if (!v.isEmpty()) nPhred = v.charAt(0) - 33;
            }
        }

        int gA = gray[1], gC = gray[2], gG = gray[3], gT = gray[4];
        char[] base = new char[256];
        int[] bandStart = new int[256];
        for (int j = 0; j < 256; ++j) { base[j] = 'N'; bandStart[j] = 0; }
        for (int j = gA; j < gC && j < 256; ++j) { base[j] = 'A'; bandStart[j] = gA; }
        for (int j = gC; j < gG && j < 256; ++j) { base[j] = 'C'; bandStart[j] = gC; }
        for (int j = gG; j < gT && j < 256; ++j) { base[j] = 'G'; bandStart[j] = gG; }
        for (int j = gT; j < gT + 63 && j < 253; ++j) { base[j] = 'T'; bandStart[j] = gT; }
        for (int j = 0; j < 256; ++j) {
            int ub = unescape(j);
            base[j] = base[ub];
            bandStart[j] = bandStart[ub];
        }

        buildInverse(qmap);
        applyDecode(qdec);

        for (int j = 0; j < 256; ++j) {
            char c = base[j];
            int val;
            if (c == 'N') val = nPhred;
            else {
                int y = unescape(j) - bandStart[j];
                if (y < 0) y = 0;
                if (y > 63) y = 63;
                val = inverse[y > 62 ? 62 : y];
            }
            if (val < 0) val = 0;
            if (val > maxPhred) val = maxPhred;
            byteLut[j] = (short) (c | ((33 + val) << 8));

            int q;
            if (c == 'N') q = nPhred;
            else {
                q = unescape(j) - bandStart[j];
                if (q < 0) q = 0;
                if (q > maxPhred) q = maxPhred;
            }
            nonetLut[j] = (short) (c | ((q & 0xFF) << 8));
        }

        char[] nb = new char[16];
        int[] ns = new int[16];
        nb[0] = 'N';
        for (int j = 1; j <= 3; ++j)   { nb[j] = 'A'; ns[j] = j - 1; }
        for (int j = 4; j <= 6; ++j)   { nb[j] = 'C'; ns[j] = j - 4; }
        for (int j = 7; j <= 9; ++j)   { nb[j] = 'G'; ns[j] = j - 7; }
        for (int j = 10; j <= 12; ++j) { nb[j] = 'T'; ns[j] = j - 10; }
        for (int j = 0; j < 256; ++j) {
            int ub = unescape(j), hi = ub >> 4, lo = ub & 0x0F, packed = 0;
            if (hi <= 12 && nb[hi] != 0) {
                int val = nb[hi] == 'N' ? nPhred : inverse[ns[hi]];
                val = Math.max(0, Math.min(maxPhred, val));
                packed |= nb[hi] | ((33 + val) << 8);
            }
            if (lo <= 12 && nb[lo] != 0) {
                int val = nb[lo] == 'N' ? nPhred : inverse[ns[lo]];
                val = Math.max(0, Math.min(maxPhred, val));
                packed |= (nb[lo] << 16) | ((33 + val) << 24);
            }
            nibLut[j] = packed;
            nibBase2[j] = (short) ((packed & 0xFF) | (((packed >>> 16) & 0xFF) << 8));
            nibQual2[j] = (short) (((packed >>> 8) & 0xFF) | (((packed >>> 24) & 0xFF) << 8));
            nibCnt[j] = (byte) (((packed & 0xFF) != 0 ? 1 : 0) + ((packed & 0xFF0000) != 0 ? 1 : 0));
        }
    }

    private void buildInverse(String qmapCsv) {
        for (int i = 0; i < 64; ++i) inverse[i] = 0;
        if (qmapCsv == null || qmapCsv.isEmpty()) return;
        int[] raw = parseIntCsv(qmapCsv, 94);
        int n = raw.length, maxSlot = 0;
        int[] lut = new int[n];
        for (int k = 0; k < n; ++k) lut[k] = Math.max(0, Math.min(63, raw[k]));
        for (int k = 0; k < n; ++k) if (lut[k] > maxSlot) maxSlot = lut[k];
        for (int k = 0; k < n; ++k) { int s = lut[k]; if (s >= 1 && s <= 63) inverse[s - 1] = k; }
        if (maxSlot >= 1 && maxSlot <= 63) {
            int rep = -1;
            for (int k = 0; k < n; ++k) if (lut[k] == maxSlot) { rep = k; break; }
            inverse[maxSlot - 1] = rep < 0 ? 0 : rep;
        }
    }

    private void applyDecode(String qdecCsv) {
        if (qdecCsv == null || qdecCsv.isEmpty()) return;
        int[] v = parseIntCsv(qdecCsv, 64);
        for (int i = 0; i < v.length; ++i) inverse[i] = Math.max(0, Math.min(93, v[i]));
    }

    private void ensure(int need) {
        if (seqBuf.length < need) {
            int m = Integer.highestOneBit(need - 1) << 1;
            seqBuf = new byte[m];
            qualBuf = new byte[m];
        }
    }

    private void readNext() throws SequenceFormatException {
        try {
            nextSequence = null;
            if (getc() < 0) return;
            if (getLine() < 0) return;
            String id = "@" + new String(line, 0, lineLen, StandardCharsets.ISO_8859_1);
            if (getLine() < 0) throw new SequenceFormatException("Truncated FASTR record after " + id);

            final byte[] raw = line;
            final int rl = lineLen;
            ensure((enc == 1 ? rl * 2 : rl) + 2);
            final byte[] sq = seqBuf, qu = qualBuf;
            int L;

            if (enc == 1) {
                L = 0;
                int nfull = rl > 0 ? rl - 1 : 0;
                for (int i = 0; i < nfull; ++i) {
                    int b = raw[i] & 0xFF;
                    short bb = nibBase2[b], qq = nibQual2[b];
                    sq[L] = (byte) bb; sq[L + 1] = (byte) (bb >> 8);
                    qu[L] = (byte) qq; qu[L + 1] = (byte) (qq >> 8);
                    L += nibCnt[b];
                }
                if (rl > 0) {
                    int v = nibLut[raw[rl - 1] & 0xFF];
                    if ((v & 0xFF) != 0)     { sq[L] = (byte) v;          qu[L] = (byte) (v >>> 8);  ++L; }
                    if ((v & 0xFF0000) != 0) { sq[L] = (byte) (v >>> 16); qu[L] = (byte) (v >>> 24); ++L; }
                }
            } else if (enc == 2) {
                int sep = -1;
                for (int i = 0; i < rl; ++i) if (raw[i] == (byte) 0xFD) { sep = i; break; }
                L = sep >= 0 ? sep : rl;
                final int mp = maxPhred;
                if (sep >= 0) {
                    int bit = sep + 1;
                    for (int k = 0; k < L; ++k) {
                        int v = nonetLut[raw[k] & 0xFF];
                        int ph = ((v >> 8) & 0xFF) + 63 * (((raw[bit + k / 7] & 0xFF) >> (k % 7)) & 1);
                        if (ph > mp) ph = mp;
                        sq[k] = (byte) v;
                        qu[k] = (byte) (33 + ph);
                    }
                } else {
                    for (int k = 0; k < L; ++k) {
                        int v = nonetLut[raw[k] & 0xFF];
                        sq[k] = (byte) v;
                        qu[k] = (byte) (33 + ((v >> 8) & 0xFF));
                    }
                }
            } else {
                L = rl;
                for (int i = 0; i < L; ++i) {
                    short v = byteLut[raw[i] & 0xFF];
                    sq[i] = (byte) v;
                    qu[i] = (byte) (v >> 8);
                }
            }

            nextSequence = new Sequence(this,
                    new String(sq, 0, L, StandardCharsets.ISO_8859_1),
                    new String(qu, 0, L, StandardCharsets.ISO_8859_1),
                    id);
        } catch (IOException e) {
            throw new SequenceFormatException(e.getMessage());
        }
    }

    @Override public String  name()         { return name; }
    @Override public boolean isColorspace() { return false; }
    @Override public File    getFile()      { return file; }
    @Override public boolean hasNext()      { return nextSequence != null; }
    public void remove() {}

    @Override public int getPercentComplete() {
        if (!hasNext() || fileSize == 0) return 100;
        try {
            return (int) Math.min(100, fis.getChannel().position() * 100 / fileSize);
        } catch (IOException e) {
            return 0;
        }
    }

    @Override public Sequence next() throws SequenceFormatException {
        Sequence s = nextSequence;
        readNext();
        return s;
    }
}
///END OF CHANGE///
