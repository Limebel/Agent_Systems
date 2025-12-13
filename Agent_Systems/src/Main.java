import java.lang.reflect.Array;
import java.util.*;
import java.io.*;
import Agents.*;

public class Main {
    public static void main(String[] args) {
        int[] w = {1,2,3};
        double[] c = {0.25,0.5,0.75};

        MFN test = new MFN.Builder()
            .m(3)
            .W(w)
            .C(c)
            .L(w)
            .R(c)
            .rho(c)
            .build();

        test.getMPs("MPs0.csv");
        test.printMPs();

    }
}