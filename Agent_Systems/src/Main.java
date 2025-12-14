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
    // reading file
        test.getMPs("MPs0.csv");
        test.printMPs();
    // cdf
        double[][] cdf = test.CDF(test.arPMFLoop());
        for (double[] row : cdf) {
            for (double value : row) {
                System.out.print(value + " ");
            }
            System.out.println();
        }
        System.out.println();
    // normalCDF
        System.out.println(test.normalCDF(1.0));
    // worstCaseNSS
        System.out.println(test.worstCaseNSS(0.01,0.05)); // 9604
        System.out.println(test.worstCaseNSS(0.05,0.05)); // 384
        System.out.println(test.worstCaseNSS(0.01,0.01)); // 16576
        System.out.println(test.worstCaseNSS(0.02,0.10)); // 1692
        System.out.println(test.worstCaseNSS(0.05,0.01)); // 664



    }
}