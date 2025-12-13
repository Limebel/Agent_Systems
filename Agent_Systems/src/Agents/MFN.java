package Agents;
import java.util.*;
import java.io.*;

public class MFN {
    private int m;
    private int[] W;
    private double[] C;
    private int[] L;
    private double[] R;
    private double[] rho;
    private double[] beta;
    private ArrayList<int[]> MPs = new ArrayList<>();

    public MFN(Builder builder) {
        this.m = builder.m;
        if (builder.W.length != m || builder.C.length != m ||
                builder.L.length != m || builder.R.length != m ||
                builder.rho.length != m) {
            throw new IllegalArgumentException("Vectors of wrong size.");
        }
        for (double value : builder.R) {
            if (value < 0 || value > 1) {
                throw new IllegalArgumentException("Invalid values in R vector.");
            }
        }
        for (double value : builder.rho) {
            if (value < 0 || value > 1) {
                throw new IllegalArgumentException("Invalid values in R vector.");
            }
        }
        this.beta = new double[m];
        for (int i = 0; i < m; i++) {
            beta[i] = 1 + (builder.rho[i] * (1 - builder.R[i])) / builder.R[i];
        }
        this.W = builder.W;
        this.C = builder.C;
        this.L = builder.L;
        this.R = builder.R;
        this.rho = builder.rho;
        this.MPs = (builder.MPs != null) ? builder.MPs : new ArrayList<>();
    }

    public static class Combinatorial{
        public static int factorial (int n){
            int result = 1;
            for (int i = 1; i <= n; i++) {
                result *= i;
            }
            return result;
        }

        public static int newtonSymbol (int a, int b) {
            if (a >= b){
                return factorial(a)/factorial(b)/factorial(a-b);
            }
            else{
                throw new IllegalArgumentException("Inwalid Newton symbol values.");
            }
        }
    }


    private double probabilityOfStateKCi_Formula1 (int k, int i){
        if (i < 0 || i >= m) throw new IllegalArgumentException("Invalid link index");
        if (k < 0 || k > W[i]) throw new IllegalArgumentException("Invalid capacity state");
        if(k > 1){
            return 1/beta[i]*Combinatorial.newtonSymbol(W[i],k)*
                    Math.pow(rho[i]*beta[i], k)*
                    Math.pow(1 - rho[i]*beta[i], W[i]-k);
        }
        else{
            return 1-1/beta[i]*
                    (1-Math.pow(1-R[i]*beta[i],W[i]));
        }
    }

    private double transmissionTimeOfd_Formula3 (int[] P, double d, double[] X){
        double capacityOfP = capacityOfPUnderX_Formula5(P, X);
        if (capacityOfP>0){
            return leadTimeOfP_Formula4(P) + Math.ceil(d / capacityOfP);
        }
        else{
            throw new IllegalArgumentException("Infinity value (Formula 3).");
        }
    }

    private int leadTimeOfP_Formula4 (int[] P){
        int sum = 0;
        for (int edge : P) {
            sum += L[edge];
        }
        return sum;
    }

    private double capacityOfPUnderX_Formula5 (int[] P, double[] X){
        double min = X[P[0]];
        for (int edge : P) {
            min = Math.min(min, X[edge]);
        }
        return min;
    }

    private double transmissionCostOfP_Formula6 (int[] P){
        double totalcost = 0;
        for (int edge : P) {
            totalcost += beta[edge];
        }
        return totalcost;
    }

    private double transmissionTimeUnderX_Formula8 (double d, double[] X, double b){
        double minTime = Double.POSITIVE_INFINITY;

        for(int[] P : MPs){
            double pathCost = transmissionCostOfP_Formula6(P);
            if (pathCost<=b){
                double time = transmissionTimeOfd_Formula3(P, d, X);
                if(time<minTime){
                    minTime = time;
                }
            }
        }
        return minTime;
    }

// void getMPs(String fileName) that reads the file with the file name = filename and creates ArrayList<int[]> MPs
    public void getMPs(String filename) {
        try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] tokens = line.trim().split(",");
                int[] nums = Arrays.stream(tokens)
                        .mapToInt(Integer::parseInt)
                        .toArray();
                MPs.add(nums);
            }
        } catch (IOException e) {
            System.err.println("Error reading file: " + e.getMessage());
        }
    }
// just for testing
    public void printMPs() {
        for (int[] row : MPs) {
            System.out.println(Arrays.toString(row));
        }
    }

// double[][] CDF(double[][] arPMF) that creates an array of values of the cumulative distribution function based on an array arPMF created by formula (1)




//  static double normalCDF(double z) that computes an approximated value of the cumulative distribution function of the standard normal distribution for n=100,
// based on the formula (https://en.wikipedia.org/wiki/Normal_distribution) ... *formula* where !! denotes the double factorial and should also be implemented

//  static double normalICDF(double u) that computes an approximated value of the quantile function (the inverse of the cumulative distribution function) of the
// standard normal distribution
// IMPORTANT! In order to implement normalICDF, invent your own algorithm that for given value u, it determines a real number x such that
// |𝑛𝑜𝑟𝑚𝑎𝑙𝐶𝐷𝐹(𝑥) − 𝑢| ≤ 10ି ଵ଴

// Based on formula (12b) (from [2] G.S. Fishman – “Monte Carlo, Concepts, Algorithms, and Applications” – Springer), implement a function finding the worst-
// case normal sample size

// Based on the inverse CDF method applied to discrete distribution or the Chen and Asau Guide Table Method coming from [3] J. E. Gentle – “Random Number
// Generation and Monte Carlo Methods” – Springer (2005), implement method double[][] randomSSV(int N, double[][]arCDF), that, for a given integer N and an
// array arCDF of values of the cumulative distribution function, generates N random system state vectors (SSVs).

    public static class Builder{
        private int m;
        private int[] W;
        private double[] C;
        private int[] L;
        private double[] R;
        private double[] rho;
        private double[] beta;
        private ArrayList<int[]> MPs;
        public Builder m(int m) {
            this.m = m;
            return this;
        }
        public Builder W(int[] W) {
            this.W = W;
            return this;
        }
        public Builder C(double[] C) {
            this.C = C;
            return this;
        }
        public Builder L(int[] L) {
            this.L = L;
            return this;
        }
        public Builder R(double[] R) {
            this.R = R;
            return this;
        }
        public Builder rho(double[] rho) {
            this.rho = rho;
            return this;
        }
        public Builder MPs(ArrayList<int[]> MPs) {
            this.MPs = MPs;
            return this;
        }
        public MFN build() {
            return new MFN(this);
        }
    }
    public static Builder builder() {
        return new Builder();
    }
}


