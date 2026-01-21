package Agents;
import jade.core.AID;
import jade.core.Agent;
import jade.domain.DFService;
import jade.domain.FIPAAgentManagement.DFAgentDescription;
import jade.domain.FIPAAgentManagement.ServiceDescription;
import jade.domain.FIPAException;
import jade.lang.acl.ACLMessage;

import java.io.IOException;
import java.util.Arrays;

public class SSVGenerator extends Agent {
    private double epsilon;
    private double delta;
    private int N;

    private MFN mfn;
    private double[][] SSVs;

    protected void setup() {
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

        N = MFN.worstCaseNSS(epsilon, delta);
        System.out.println("The minimum number of iterations is equal to " + N);

        SSVGeneratorGui gui = new SSVGeneratorGui(this);
        gui.setVisible(true);
    }

    public void onSendData(
            int[] W,
            double[] C,
            int[] L,
            double[] R,
            double[] rho,
            String mpFile
    ){
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

        double[][] arPMF = mfn.arPMFLoop();
        double[][] arCDF = mfn.CDF(arPMF);
        SSVs = mfn.randomSSV(N, arCDF);
        System.out.println(N + " random SSVs have been generated!");

        for (int i = 0; i < 5; i++) {
            System.out.println(Arrays.toString(SSVs[i]));
        }

        AID ttAgent = findTTAgent();
        if (ttAgent == null) {
            System.out.println("TT agent not found!");
            doDelete();
            return;
        }
        sendDataToTT(ttAgent, W, C, L, R, rho, SSVs, mpFile);
    }

    private AID findTTAgent() {
        DFAgentDescription template = new DFAgentDescription();
        ServiceDescription sd = new ServiceDescription();
        sd.setType("transmission-time");
        template.addServices(sd);

        try {
            DFAgentDescription[] result =
                    DFService.search(this, template);

            if (result.length > 0) {
                return result[0].getName();
            }
        } catch (FIPAException e) {
            e.printStackTrace();
        }

        return null;
    }

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
}
