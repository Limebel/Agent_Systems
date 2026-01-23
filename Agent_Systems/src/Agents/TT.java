package Agents;

import jade.core.AID;
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
    private double T;   // max transmission time we want

    @Override
    protected void setup() {
        readArguments(); //Initialization method
        registerService(); //Registering service
        waitForSSVMessage(); //Waiting to be chosen for the task
    }

    //Initialization method, just loading and  checking arguments correctness
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

    //Informing on some kind of service registering page(dfd) that I want to do specific service
    private void registerService() {
        DFAgentDescription dfd = new DFAgentDescription();
        dfd.setName(getAID());

        ServiceDescription sd = new ServiceDescription();
        sd.setType("transmission-time"); //Type of service
        sd.setName("TT-agent");

        dfd.addServices(sd);
        System.out.println("TT service registered");
        try {
            DFService.register(this, dfd);
        } catch (FIPAException e) {
            e.printStackTrace();
        }
    }

    //Waiting until receiving information with task
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

    //Doing the task
    private void handleMessage(ACLMessage msg) {
        try {
            //Unpacking message
            Object[] data = (Object[]) msg.getContentObject();

            int[] W = (int[]) data[0];
            double[] C = (double[]) data[1];
            int[] L = (int[]) data[2];
            double[] R = (double[]) data[3];
            double[] rho = (double[]) data[4];
            double[][] SSVs = (double[][]) data[5];
            String mpFilePath = (String) data[6];

            System.out.println("Received SSV data from SSVGenerator");

            writeSSVsToCSV(SSVs); //Writing SSVs to csv file (no needed for the task, required by instruction)
            double reliability = computeReliability(W, C, L, R, rho, SSVs, mpFilePath); //Calculating reliability

            sendResult(msg.getSender(), reliability); //Sending results back

        } catch (Exception e) {
            e.printStackTrace();
        }

        doDelete(); //Deletion after oing the task
    }

    //Writing SSVs into csv file
    private void writeSSVsToCSV(double[][] SSVs){
        File csvOutputFile = new File("SSV.csv");
        try(PrintWriter pw = new PrintWriter(csvOutputFile)){
            for(double[] ssv : SSVs){
                for (int i = 0; i < ssv.length; i++) {
                    pw.print((int) ssv[i]);
                    if (i < ssv.length - 1) pw.print(",");
                }
                pw.println();
            }
        } catch(Exception e){
            e.printStackTrace();
        }
    }

    //Computing reliability (main task)
    private double computeReliability(int[] W, double[] C, int[] L, double[] R, double[]rho, double[][] SSVs, String mpFilePath){
        //Creating mfn object based on information
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

        //Reliability is a probability that we can send amount of flow (d) through system (mfn) in shorter time than (T)
        //We calculate it by just testing many different states of a system (SSVs)
        for(double[] X : SSVs){
            if(mfn.transmissionTimeUnderX_Formula8(d, X) <= T){
                success++;
            }
        }
        return (double) success/SSVs.length;
    }

    //Sending the results to SSVGenerator agent
    private void sendResult(AID receiver, double reliability) {
        ACLMessage reply = new ACLMessage(ACLMessage.INFORM);
        reply.addReceiver(receiver);
        reply.setContent("Estimated reliability = " + reliability);
        send(reply);
    }
}
