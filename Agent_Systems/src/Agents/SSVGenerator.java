package Agents;
import jade.core.AID;
import jade.core.Agent;
import jade.core.behaviours.CyclicBehaviour;
import jade.domain.DFService;
import jade.domain.FIPAAgentManagement.DFAgentDescription;
import jade.domain.FIPAAgentManagement.ServiceDescription;
import jade.domain.FIPAException;
import jade.lang.acl.ACLMessage;

import java.io.IOException;
import java.util.Arrays;

public class SSVGenerator extends Agent {
    private double epsilon; //accuracy we want
    private double delta; //confidence we want
    private int N; //number of SSVs we want (calculate using worst case NSS)

    private MFN mfn;
    private double[][] SSVs;

    protected void setup() {
        //Passing and checking correctness of epsilon and delta arguments
        Object[] args = getArguments();
        if (args == null || args.length != 2) {
            doDelete();
            return;
        }

        epsilon = Double.parseDouble(args[0].toString());
        delta   = Double.parseDouble(args[1].toString());

        if (epsilon <= 0 || epsilon >= 1 || delta <= 0 || delta >= 1) {
            System.out.println("ε and δ must be in (0,1)");
            doDelete();
            return;
        }

        //Calculating number of SSVs needed
        N = MFN.worstCaseNSS(epsilon, delta);
        System.out.println("The minimum number of iterations is equal to " + N);

        //Activating SSVGeneratorGUI to pass other arguments about System of links
        SSVGeneratorGui gui = new SSVGeneratorGui(this);
        gui.setVisible(true);
    }

    //Method that executes after data about system(mfn) is correctly send by user
    public void onSendData(
            int[] W,
            double[] C,
            int[] L,
            double[] R,
            double[] rho,
            String mpFile
    ){
        //Creating mfn object
        mfn = MFN.builder()
                .m(W.length)
                .W(W)
                .C(C)
                .L(L)
                .R(R)
                .rho(rho)
                .build();
        mfn.getMPs(mpFile);
        System.out.println("It has been created the MFN with the following parameters:");
        System.out.println("W=" + Arrays.toString(W));
        System.out.println("C=" + Arrays.toString(C));
        System.out.println("L=" + Arrays.toString(L));
        System.out.println("R=" + Arrays.toString(R));
        System.out.println("rho=" + Arrays.toString(rho));

        //Calculating probabilities used for generating SSVs
        double[][] arPMF = mfn.arPMFLoop();
        double[][] arCDF = mfn.CDF(arPMF);

        //Generating random SSVs
        SSVs = mfn.randomSSV(N, arCDF);
        System.out.println(N + " random SSVs have been generated!");

        //for printing SSVs
        /*for (int i = 0; i < 5; i++) {
            System.out.println(Arrays.toString(SSVs[i]));
        }*/

        AID ttAgent = findTTAgent(); //Trying to find agent to calculate reliability
        if (ttAgent == null) {
            System.out.println("TT agent not found!");
            doDelete();
            return;
        }
        sendDataToTT(ttAgent, W, C, L, R, rho, SSVs, mpFile); //sending data needed for the task to agent
        waitForTTReply(); //waiting for a reply with task results

    }

    //Method for finding agent to calculate reliability
    private AID findTTAgent() {
        DFAgentDescription template = new DFAgentDescription();
        ServiceDescription sd = new ServiceDescription();
        sd.setType("transmission-time"); //SSVGenerator indicates that it looks for TT agent
        template.addServices(sd);
        try {
            //Collecting all agents that requested doing the task
            DFAgentDescription[] result =
                    DFService.search(this, template);// Searching through service registering board

            //Choosing the first agent to do the task
            if (result.length > 0) {
                System.out.println("Agent" + result[0].getName() + "found");
                return result[0].getName();
            }
        } catch (FIPAException e) {
            e.printStackTrace();
        }

        return null;
    }

    //Function for sending data needed for task to agent
    private void sendDataToTT(AID ttAgent, int[] W, double[] C, int[] L, double[] R, double[]rho, double[][] SSVs, String mpFilePath) {
        try {
            ACLMessage msg = new ACLMessage(ACLMessage.INFORM);
            msg.addReceiver(ttAgent);
            msg.setContentObject(new Object[] {W, C, L, R, rho, SSVs, mpFilePath});
            send(msg);

            System.out.println("Data sent to TT agent");

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    //Waiting for a reply. After getting message printing it and deleting itself
    private void waitForTTReply() {
        addBehaviour(new CyclicBehaviour() {
            @Override
            public void action() {
                ACLMessage msg = receive();
                if (msg != null) {
                    System.out.println(msg.getContent());
                    doDelete();
                } else {
                    block();
                }
            }
        });
    }
}
