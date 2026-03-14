import java.io.*;
import java.util.*;

public class FASTRDebug {

    static final int SENTINEL = 0xFF;
    static final int PHRED_MAX = 93;
    static final int PHRED_OFFSET = 33;

    static char[] BASE_FOR_VALUE  = new char[256];
    static int[]  LOWER_FOR_VALUE = new int[256];
    static int[]  INVERSE_SCALE   = new int[64];

    public static void main(String[] args) throws Exception {
        if (args.length < 1) { System.err.println("Usage: java FASTRDebug <file.fastr>"); return; }

        File file = new File(args[0]);
        byte[] data = new FileInputStream(file).readAllBytes();

        // parse header
        int gN=0, gA=3, gG=66, gC=129, gT=192;
        int phredMax = PHRED_MAX;
        int cursor = 0;
        while (cursor < data.length) {
            int lineEnd = indexOf(data, (byte)'\n', cursor);
            if (lineEnd == -1) lineEnd = data.length;
            String line = new String(data, cursor, lineEnd - cursor, "UTF-8").trim();
            cursor = lineEnd + 1;
            if (line.startsWith("@")) { cursor = lineEnd - line.length(); break; }
            if (line.startsWith("#GRAY_VALS=")) {
                String gv = line.substring(11).replaceAll("[\\[\\]\\s]", "");
                String[] p = gv.split(",");
                if (p.length >= 5) { gN=Integer.parseInt(p[0]); gA=Integer.parseInt(p[1]); gG=Integer.parseInt(p[2]); gC=Integer.parseInt(p[3]); gT=Integer.parseInt(p[4]); }
                System.out.println("GRAY_VALS: N="+gN+" A="+gA+" G="+gG+" C="+gC+" T="+gT);
            }
            if (line.startsWith("#PHRED-ALPHABET=")) {
                String pa = line.substring(16).trim();
                if (pa.startsWith("PHRED_")) phredMax = Integer.parseInt(pa.substring(6)) - 1;
                System.out.println("PHRED_MAX: " + phredMax);
            }
            if (line.startsWith("#QUAL_SCALE=")) {
                System.out.println("QUAL_SCALE: " + line.substring(12));
            }
            if (line.startsWith("#MODE=")) {
                System.out.println("MODE: " + line.substring(6));
            }
        }

        buildLookupTables(gN, gA, gG, gC, gT);
        buildInverseScaleTable(phredMax);

        System.out.println("\n--- INVERSE_SCALE table (scaled -> phred) ---");
        for (int s = 0; s < 64; s++) System.out.printf("s=%2d -> phred=%d%n", s, INVERSE_SCALE[s]);

        // find first record
        int headerStart = cursor;
        while (headerStart < data.length && data[headerStart] != '@') headerStart++;
        int xffPos = indexOf(data, (byte)0xFF, headerStart);
        if (xffPos == -1) { System.err.println("No sentinel found"); return; }

        String header = new String(data, headerStart+1, xffPos-headerStart-1, "UTF-8").trim();
        System.out.println("\n--- First record header: " + header);

        int seqStart = xffPos + 1;
        if (seqStart < data.length && data[seqStart] == '\n') seqStart++;

        // find body end (next \n@\xff pattern)
        int seqEnd = findBodyEnd(data, seqStart);
        int bodyLen = seqEnd - seqStart;
        System.out.println("Body length: " + bodyLen + " bytes");

        // print first 20 bytes decoded
        System.out.println("\n--- First 20 bytes ---");
        System.out.printf("%-6s %-6s %-6s %-8s %-8s %-8s%n", "idx", "raw", "base", "lower", "scaled", "phred");
        for (int i = 0; i < Math.min(20, bodyLen); i++) {
            int v      = data[seqStart + i] & 0xFF;
            char base  = BASE_FOR_VALUE[v];
            int lower  = LOWER_FOR_VALUE[v];
            int scaled = Math.max(0, Math.min(63, v - lower));
            int phred  = INVERSE_SCALE[scaled];
            System.out.printf("%-6d %-6d %-6c %-8d %-8d %-8d%n", i, v, base, lower, scaled, phred);
        }
    }

    static int findBodyEnd(byte[] data, int seqStart) {
        int candidateEnd = seqStart;
        while (true) {
            int nextAt = indexOfSeq(data, new byte[]{'\n','@'}, candidateEnd);
            if (nextAt == -1) {
                int end = data.length;
                while (end > seqStart && (data[end-1]=='\n'||data[end-1]=='\r'||data[end-1]==' ')) end--;
                return end;
            }
            int nextXff     = indexOf(data, (byte)0xFF, nextAt+1);
            int followingAt = indexOfSeq(data, new byte[]{'\n','@'}, nextAt+2);
            if (followingAt != -1 && (nextXff == -1 || followingAt < nextXff)) {
                candidateEnd = followingAt + 1;
                continue;
            }
            return nextAt;
        }
    }

    static void buildLookupTables(int gN, int gA, int gG, int gC, int gT) {
        for (int v=0;v<256;v++){BASE_FOR_VALUE[v]='N';LOWER_FOR_VALUE[v]=0;}
        for (int v=gN;v<gA&&v<255;v++){BASE_FOR_VALUE[v]='N';LOWER_FOR_VALUE[v]=gN;}
        for (int v=gA;v<gG&&v<255;v++){BASE_FOR_VALUE[v]='A';LOWER_FOR_VALUE[v]=gA;}
        for (int v=gG;v<gC&&v<255;v++){BASE_FOR_VALUE[v]='G';LOWER_FOR_VALUE[v]=gG;}
        for (int v=gC;v<gT&&v<255;v++){BASE_FOR_VALUE[v]='C';LOWER_FOR_VALUE[v]=gC;}
        for (int v=gT;v<255;v++)      {BASE_FOR_VALUE[v]='T';LOWER_FOR_VALUE[v]=gT;}
    }

    static void buildInverseScaleTable(int phredMax) {
        double ln63 = Math.log(63.0);
        int[] scaledIntForQ = new int[phredMax + 1];
        for (int q = 0; q <= phredMax; q++) {
            double s = 1.0 + 62.0 * (Math.log(q + 1.0) / ln63);
            scaledIntForQ[q] = (int) s;
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
            for (int i=0;i<valid[0];i++) inverse[i]=inverse[valid[0]];
            for (int i=0;i<nValid-1;i++) {
                int s=valid[i],e=valid[i+1],sv=inverse[s],ev=inverse[e],gap=e-s;
                for (int j=1;j<gap;j++) inverse[s+j]=(int)Math.round(sv+(double)(ev-sv)*j/gap);
            }
            for (int i=valid[nValid-1]+1;i<64;i++) inverse[i]=inverse[valid[nValid-1]];
        }
        for (int s=0;s<64;s++) {
            int v=inverse[s];
            INVERSE_SCALE[s]=Math.max(0,Math.min(phredMax,v==-1?0:v));
        }
    }

    static int indexOf(byte[] data, byte target, int from) {
        for (int i=from;i<data.length;i++) if(data[i]==target) return i;
        return -1;
    }
    static int indexOfSeq(byte[] data, byte[] seq, int from) {
        outer: for (int i=from;i<=data.length-seq.length;i++) {
            for (int j=0;j<seq.length;j++) if(data[i+j]!=seq[j]) continue outer;
            return i;
        }
        return -1;
    }
}
