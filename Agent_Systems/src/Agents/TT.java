package Agents;

import jade.core.Agent;
import jade.core.behaviours.CyclicBehaviour;
import jade.domain.DFService;
import jade.domain.FIPAAgentManagement.DFAgentDescription;
import jade.domain.FIPAAgentManagement.ServiceDescription;
import jade.domain.FIPAException;
import jade.lang.acl.ACLMessage;

import java.io.File;
import java.io.PrintWriter;

public class TT extends Agent {
    private double d;   // demanded flow
    private double T;   // max transmission time

    @Override
    protected void setup() {
        readArguments();
        registerService();
        waitForSSVMessage();
    }

    private void readArguments() {
        Object[] args = getArguments();

        if (args == null || args.length != 2) {
            System.out.println("TT requires two arguments: d and T");
            doDelete();
            return;
        }

        try {
            d = Double.parseDouble(args[0].toString());
            T = Double.parseDouble(args[1].toString());
        } catch (NumberFormatException e) {
            System.out.println("Invalid arguments for TT agent");
            doDelete();
        }
    }

    private void registerService() {
        DFAgentDescription dfd = new DFAgentDescription();
        dfd.setName(getAID());

        ServiceDescription sd = new ServiceDescription();
        sd.setType("transmission-time");
        sd.setName("TT-agent");

        dfd.addServices(sd);
        System.out.println("TT service registered");
        try {
            DFService.register(this, dfd);
        } catch (FIPAException e) {
            e.printStackTrace();
        }
    }

    private void waitForSSVMessage() {
        addBehaviour(new CyclicBehaviour() {
            @Override
            public void action() {
                ACLMessage msg = receive();
                if (msg != null) {
                    handleMessage(msg);
                    removeBehaviour(this);
                } else {
                    block();
                }
            }
        });
    }

    private void handleMessage(ACLMessage msg) {
        try {
            Object[] data = (Object[]) msg.getContentObject();

            int[] W = (int[]) data[0];
            double[] C = (double[]) data[1];
            int[] L = (int[]) data[2];
            double[] R = (double[]) data[3];
            double[] rho = (double[]) data[4];
            double[][] SSVs = (double[][]) data[5];
            String mpFilePath = (String) data[6];

            System.out.println("Received SSV data from SSVGenerator");

            writeSSVsToCSV(SSVs);
            double reliability = computeReliability(W, C, L, R, rho, SSVs, mpFilePath);

            //sendResult(msg.getSender(), reliability);

        } catch (Exception e) {
            e.printStackTrace();
        }

        doDelete();
    }

    private void writeSSVsToCSV(double[][] SSVs){
        File csvOutputFile = new File("SSV.csv");
        try(PrintWriter pw = new PrintWriter(csvOutputFile)){
            for(double[] ssv : SSVs){
                for(double state : ssv){
                    pw.print((int)state);
                    pw.print(",");
                }
                pw.println();
            }
        } catch(Exception e){
            e.printStackTrace();
        }
    }

    private double computeReliability(int[] W, double[] C, int[] L, double[] R, double[]rho, double[][] SSVs, String mpFilePath){
        MFN mfn = MFN.builder()
                .m(W.length)
                .W(W)
                .C(C)
                .L(L)
                .R(R)
                .rho(rho)
                .build();

        mfn.getMPs(mpFilePath);

        int success = 0;

        for(double[] X : SSVs){
            if(mfn.transmissionTimeUnderX_Formula8(d, X) < T){
                success++;
            }
        }
        return (double) success/SSVs.length;
    }
}
