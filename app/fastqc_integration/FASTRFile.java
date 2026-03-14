package uk.ac.babraham.FastQC.Sequence;

import java.io.*;
import java.util.*;
import javax.script.*;

import uk.ac.babraham.FastQC.FastQCConfig;

public class FASTRFile implements SequenceFile {

    private static final int SENTINEL = 0xFF;

    private final File file;
    private final long fileSize;
    private final String name;
    private final FileInputStream fis;

    private Sequence nextSequence = null;
    private long seqIndex = 1;

    private int mode = 2;
    private int phredOffset = 33;
    private int phredMax = 93;
    private String sraAccession = null;
    private String structureTemplate = null;
    private String qualScaleFormula = "x";

    private String[] headersTable = null;

    private final char[] BASE_FOR_VALUE  = new char[256];
    private final int[]  LOWER_FOR_VALUE = new int[256];
    private final int[]  INVERSE_SCALE   = new int[64];

    private byte[] data;
    private int cursor = 0;

    protected FASTRFile(FastQCConfig config, File file) throws SequenceFormatException, IOException {
        this.file     = file;
        this.fileSize = file.length();
        this.name     = file.getName();

        fis = new FileInputStream(file);
        data = fis.readAllBytes();

        int gN = 0, gA = 3, gG = 66, gC = 129, gT = 192;
        Map<String, String> meta = new HashMap<>();

        int pos = 0;
        while (pos < data.length) {
            int lineEnd = indexOf(data, (byte)'\n', pos);
            if (lineEnd == -1) lineEnd = data.length;
            String line = new String(data, pos, lineEnd - pos, "UTF-8").trim();
            pos = lineEnd + 1;
            if (line.startsWith("@")) { cursor = lineEnd - line.length(); break; }
            if (line.startsWith("#") && line.contains("=")) {
                String[] kv = line.substring(1).split("=", 2);
                meta.put(kv[0].trim().toUpperCase(), kv[1].trim());
            }
        }

        if (meta.containsKey("MODE")) {
            try { mode = Integer.parseInt(meta.get("MODE")); } catch (NumberFormatException ignored) {}
        }
        if (meta.containsKey("ACCESSION")) {
            String acc = meta.get("ACCESSION");
            if (!acc.isEmpty()) sraAccession = acc;
        }
        if (meta.containsKey("PHRED-ALPHABET")) {
            String pa = meta.get("PHRED-ALPHABET");
            if (pa.startsWith("PHRED_")) {
                try { phredMax = Integer.parseInt(pa.substring(6)) - 1; } catch (NumberFormatException ignored) {}
            }
        }
        if (meta.containsKey("STRUCTURE")) {
            structureTemplate = meta.get("STRUCTURE");
        }
        if (meta.containsKey("GRAY_VALS")) {
            String gv = meta.get("GRAY_VALS").replaceAll("[\\[\\]\\s]", "");
            String[] parts = gv.split(",");
            if (parts.length >= 5) {
                try {
                    gN = Integer.parseInt(parts[0]);
                    gA = Integer.parseInt(parts[1]);
                    gG = Integer.parseInt(parts[2]);
                    gC = Integer.parseInt(parts[3]);
                    gT = Integer.parseInt(parts[4]);
                } catch (NumberFormatException ignored) {}
            }
        }
        if (meta.containsKey("QUAL_SCALE") && !meta.get("QUAL_SCALE").isEmpty()) {
            qualScaleFormula = convertFormula(meta.get("QUAL_SCALE"));
        }

        buildLookupTables(gN, gA, gG, gC, gT);
        buildInverseScaleTable();
        if (mode == 3) loadSidecarHeaders();
        readNext();
    }

    private String convertFormula(String formula) {

        String f = formula.trim().replaceAll("^\\s*f\\s*\\(\\s*x\\s*\\)\\s*=\\s*", "");
        f = f.replace("^", "**");
        f = f.replaceAll("\\bln\\s*\\(", "Math.log(");
        f = f.replaceAll("\\blog10\\s*\\(", "Math.log10(");
        f = f.replaceAll("\\blog\\s*\\(", "Math.log(");
        f = f.replaceAll("\\bexp\\s*\\(", "Math.exp(");
        f = f.replaceAll("\\bsqrt\\s*\\(", "Math.sqrt(");
        f = f.replaceAll("\\babs\\s*\\(", "Math.abs(");
        return f;
    }

    @Override public String  name()         { return name; }
    @Override public boolean isColorspace() { return false; }
    @Override public File    getFile()      { return file; }
    @Override public boolean hasNext()      { return nextSequence != null; }

    @Override
    public Sequence next() throws SequenceFormatException {
        Sequence s = nextSequence;
        readNext();
        return s;
    }

    @Override
    public int getPercentComplete() {
        if (!hasNext()) return 100;
        return (int)((double) cursor / data.length * 100);
    }

    private void readNext() throws SequenceFormatException {
        try {
            if      (mode == 0)              readNextMode0();
            else if (mode == 1 || mode == 2) readNextMode1or2();
            else if (mode == 3)              readNextMode3();
            else throw new SequenceFormatException("Unknown FASTR mode: " + mode);
        } catch (IOException ioe) {
            throw new SequenceFormatException(ioe.getMessage());
        }
    }

    private void readNextMode0() throws IOException, SequenceFormatException {
        while (cursor < data.length && data[cursor] != '@') cursor++;
        if (cursor >= data.length) { nextSequence = null; return; }

        int line1End = indexOf(data, (byte)'\n', cursor);
        if (line1End == -1) { nextSequence = null; return; }
        String miniHeader = new String(data, cursor + 1, line1End - cursor - 1, "UTF-8").trim();
        cursor = line1End + 1;

        int line2End = indexOf(data, (byte)'\n', cursor);
        if (line2End == -1) throw new SequenceFormatException("Mode 0: truncated seq");
        String seq = new String(data, cursor, line2End - cursor, "UTF-8").trim();
        cursor = line2End + 1;

        int line3End = indexOf(data, (byte)'\n', cursor);
        if (line3End == -1) throw new SequenceFormatException("Mode 0: truncated mid");
        cursor = line3End + 1;

        int qualEnd = cursor + seq.length();
        if (qualEnd > data.length) throw new SequenceFormatException("Mode 0: truncated quality");
        String quality = new String(data, cursor, seq.length(), "UTF-8");
        cursor = qualEnd;
        if (cursor < data.length && data[cursor] == '\n') cursor++;

        String id = reconstructHeader(miniHeader);
        nextSequence = new Sequence(this, seq.toUpperCase(), quality, id);
        seqIndex++;
    }

    private void readNextMode1or2() throws IOException, SequenceFormatException {
        int headerStart = findValidatedHeader(cursor);
        if (headerStart == -1) { nextSequence = null; return; }

        int xffPos = indexOf(data, (byte)0xFF, headerStart);
        if (xffPos == -1) { nextSequence = null; return; }

        String headerContent = new String(data, headerStart + 1, xffPos - headerStart - 1, "UTF-8").trim();

        int seqStart = xffPos + 1;
        if (seqStart < data.length && data[seqStart] == '\n') seqStart++;

        int seqEnd = findBodyEnd(seqStart);

        byte[] body = new byte[seqEnd - seqStart];
        System.arraycopy(data, seqStart, body, 0, body.length);
        cursor = seqEnd;

        String id      = reconstructHeader(headerContent);
        String seq     = decodeSequence(body);
        String quality = decodeQuality(body);

        nextSequence = new Sequence(this, seq, quality, id);
        seqIndex++;
    }

    private void readNextMode3() throws IOException, SequenceFormatException {
        int xffPos = indexOf(data, (byte)0xFF, cursor);
        if (xffPos == -1) { nextSequence = null; return; }
        int seqStart = xffPos + 1;

        int nextXff = indexOf(data, (byte)0xFF, seqStart);
        int seqEnd = (nextXff == -1) ? data.length : nextXff;
        while (seqEnd > seqStart && data[seqEnd - 1] == '\n') seqEnd--;

        if (seqEnd <= seqStart) { cursor = seqEnd; readNextMode3(); return; }

        byte[] body = new byte[seqEnd - seqStart];
        System.arraycopy(data, seqStart, body, 0, body.length);
        cursor = seqEnd;

        String id;
        int idx = (int)(seqIndex - 1);
        if (headersTable != null && idx < headersTable.length) {
            String h = headersTable[idx];
            id = h.startsWith("@") ? h : "@" + h;
        } else if (sraAccession != null) {
            id = "@" + sraAccession + "." + seqIndex;
        } else {
            id = "@seq" + seqIndex;
        }

        nextSequence = new Sequence(this, decodeSequence(body), decodeQuality(body), id);
        seqIndex++;
    }

    private void loadSidecarHeaders() {
        String baseName = file.getName();
        int dot = baseName.lastIndexOf('.');
        String stem = (dot >= 0) ? baseName.substring(0, dot) : baseName;
        File dir = file.getParentFile();
        if (dir == null) dir = new File(".");

        File sidecar = new File(dir, stem + "_headers.txt");
        if (!sidecar.exists()) return;

        try {
            List<String> lines = new ArrayList<>();
            BufferedReader br = new BufferedReader(new FileReader(sidecar));
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (!line.isEmpty()) lines.add(line);
            }
            br.close();
            headersTable = lines.toArray(new String[0]);
        } catch (IOException ignored) {}
    }

    private int findValidatedHeader(int fromPos) {
        int candidate;
        if (fromPos < data.length && data[fromPos] == '@') {
            candidate = fromPos;
        } else {
            int na = indexOfSeq(data, new byte[]{'\n', '@'}, fromPos);
            if (na == -1) return -1;
            candidate = na + 1;
        }
        while (candidate >= 0 && candidate < data.length) {
            int nextXff = indexOf(data, (byte)0xFF, candidate);
            int nextAt  = indexOfSeq(data, new byte[]{'\n', '@'}, candidate + 1);
            if (nextAt != -1 && (nextXff == -1 || nextAt < nextXff)) {
                candidate = nextAt + 1;
                continue;
            }
            if (nextXff != -1) return candidate;
            return -1;
        }
        return -1;
    }

    private int findBodyEnd(int seqStart) {
        int candidateEnd = seqStart;
        while (true) {
            int nextAt = indexOfSeq(data, new byte[]{'\n', '@'}, candidateEnd);
            if (nextAt == -1) {
                int end = data.length;
                while (end > seqStart && (data[end-1] == '\n' || data[end-1] == '\r' || data[end-1] == ' ')) end--;
                return end;
            }
            int nextXff     = indexOf(data, (byte)0xFF, nextAt + 1);
            int followingAt = indexOfSeq(data, new byte[]{'\n', '@'}, nextAt + 2);
            if (followingAt != -1 && (nextXff == -1 || followingAt < nextXff)) {
                candidateEnd = followingAt + 1;
                continue;
            }
            return nextAt;
        }
    }

    private String decodeSequence(byte[] body) {
        char[] bases = new char[body.length];
        for (int i = 0; i < body.length; i++)
            bases[i] = BASE_FOR_VALUE[body[i] & 0xFF];
        return new String(bases);
    }

    private String decodeQuality(byte[] body) {
        char[] qual = new char[body.length];
        for (int i = 0; i < body.length; i++) {
            int v      = body[i] & 0xFF;
            int scaled = Math.max(0, Math.min(63, v - LOWER_FOR_VALUE[v]));
            int phred  = Math.max(0, Math.min(phredMax, INVERSE_SCALE[scaled]));
            qual[i]    = (char)(phred + phredOffset);
        }
        return new String(qual);
    }

    private String reconstructHeader(String miniContent) {
        if (structureTemplate != null && !structureTemplate.isEmpty()) {
            String[] fields = miniContent.split("[:\\s/]+");
            String result = structureTemplate;
            for (int i = 0; i < fields.length; i++)
                result = result.replace("{REPEATING_" + (i + 1) + "}", fields[i]);
            result = result.replaceAll("\\{REPEATING_\\d+\\}", "");
            return "@" + result.trim();
        }
        if (sraAccession != null && !miniContent.startsWith(sraAccession))
            return "@" + sraAccession + "." + seqIndex + " " + miniContent;
        return "@" + miniContent;
    }

    private void buildLookupTables(int gN, int gA, int gG, int gC, int gT) {
        for (int v = 0; v < 256; v++) { BASE_FOR_VALUE[v] = 'N'; LOWER_FOR_VALUE[v] = 0; }
        for (int v = gN; v < gA && v < 255; v++) { BASE_FOR_VALUE[v] = 'N'; LOWER_FOR_VALUE[v] = gN; }
        for (int v = gA; v < gG && v < 255; v++) { BASE_FOR_VALUE[v] = 'A'; LOWER_FOR_VALUE[v] = gA; }
        for (int v = gG; v < gC && v < 255; v++) { BASE_FOR_VALUE[v] = 'G'; LOWER_FOR_VALUE[v] = gG; }
        for (int v = gC; v < gT && v < 255; v++) { BASE_FOR_VALUE[v] = 'C'; LOWER_FOR_VALUE[v] = gC; }
        for (int v = gT; v < 255;      v++) { BASE_FOR_VALUE[v] = 'T'; LOWER_FOR_VALUE[v] = gT; }
    }

    private void buildInverseScaleTable() throws SequenceFormatException {
        ScriptEngine engine = null;
        try {
            ScriptEngineManager manager = new ScriptEngineManager();
            engine = manager.getEngineByName("JavaScript");
        } catch (Exception ignored) {}

        int[] scaledIntForQ = new int[phredMax + 1];
        for (int q = 0; q <= phredMax; q++) {
            double result = evalFormula(engine, qualScaleFormula, q);

            if (Double.isNaN(result)) result = 1.0;
            else if (Double.isInfinite(result)) result = result > 0 ? 63.0 : 1.0;
            scaledIntForQ[q] = (int) result;
        }

        int[] inverse = new int[64];
        Arrays.fill(inverse, -1);
        for (int q = 0; q <= phredMax; q++) {
            int s = scaledIntForQ[q];
            if (s >= 0 && s <= 63) inverse[s] = q;
        }

        int[] valid = new int[64]; int nValid = 0;
        for (int i = 0; i < 64; i++) if (inverse[i] != -1) valid[nValid++] = i;

        if (nValid > 0) {
            for (int i = 0; i < valid[0]; i++) inverse[i] = inverse[valid[0]];
            for (int i = 0; i < nValid - 1; i++) {
                int s = valid[i], e = valid[i+1];
                int sv = inverse[s], ev = inverse[e], gap = e - s;
                for (int j = 1; j < gap; j++)
                    inverse[s+j] = (int) Math.round(sv + (double)(ev-sv)*j/gap);
            }
            for (int i = valid[nValid-1]+1; i < 64; i++) inverse[i] = inverse[valid[nValid-1]];
        } else {
            Arrays.fill(inverse, 0);
        }

        for (int s = 0; s < 64; s++) {
            int v = inverse[s];
            INVERSE_SCALE[s] = Math.max(0, Math.min(phredMax, v == -1 ? 0 : v));
        }
    }

    private double evalFormula(ScriptEngine engine, String formula, double x) {
        if (engine == null) {

            if (formula.trim().equals("x")) return x;
            double ln63 = Math.log(63.0);
            return 1.0 + 62.0 * (Math.log(x + 1.0) / ln63);
        }
        try {
            engine.put("x", x);

            engine.eval("var ln=Math.log, log=Math.log, log10=Math.log10, exp=Math.exp, sqrt=Math.sqrt, abs=Math.abs;");
            Object result = engine.eval(formula);
            return ((Number) result).doubleValue();
        } catch (Exception e) {

            return x;
        }
    }

    private static int indexOf(byte[] data, byte target, int from) {
        for (int i = from; i < data.length; i++) if (data[i] == target) return i;
        return -1;
    }

    private static int indexOfSeq(byte[] data, byte[] seq, int from) {
        outer:
        for (int i = from; i <= data.length - seq.length; i++) {
            for (int j = 0; j < seq.length; j++) if (data[i+j] != seq[j]) continue outer;
            return i;
        }
        return -1;
    }

    public void remove() {}
}
