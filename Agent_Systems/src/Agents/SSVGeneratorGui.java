package Agents;

import javax.swing.*;
import java.awt.*;
import java.io.File;

public class SSVGeneratorGui extends JFrame {
    private SSVGenerator agent;

    private JTextField wField;
    private JTextField cField;
    private JTextField lField;
    private JTextField rField;
    private JTextField rhoField;

    private JTextField fileField;
    private JButton browseButton;
    private JButton sendButton;

    public SSVGeneratorGui(SSVGenerator agent) {
        this.agent = agent;
        initComponents();
    }

    private void initComponents() {
        //Setting window for passing arguments
        setTitle("SSV Generator");
        setSize(600, 400);
        setLocationRelativeTo(null);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        setLayout(new GridLayout(7, 2, 5, 5));

        // Creating fields with default values of System for example 1
        wField   = new JTextField("4,3,2,3,2");
        cField   = new JTextField("10,15,25,15,20");
        lField   = new JTextField("5,7,6,5,8");
        rField   = new JTextField("0.7,0.65,0.67,0.71,0.75");
        rhoField = new JTextField("0.1,0.3,0.5,0.7,0.9");

        //values for example. Not working, because of exceeding int limit (47!)
        /*wField   = new JTextField("43, 32, 41, 19, 15, 12, 29, 27, 20, 11, 45, 25, 25, 37, 46, 17, 26," +
                " 16, 28, 27, 33, 30, 47, 23, 31, 10, 20, 28, 31, 36, 12, 44, 26, 41," +
                " 38, 16, 21, 28, 35, 37, 42");
        cField   = new JTextField("6.25, 14.65, 9.55, 13.95, 12.0 , 6.5 , 6.89, 9.36, 13.72," +
                " 7.9 , 6.66, 8.14, 13.84, 10.53, 7.28, 8.53, 8.92, 10.22," +
                " 5.39, 7.09, 12.2 , 11.32, 10.15, 7.7 , 11.25, 9.2 , 13.91," +
                " 8.2 , 13.45, 8.54, 13.34, 5.17, 7.7 , 7.44, 8.74, 14.01," +
                " 5.53, 10.11, 6.39, 12.59, 6.03");
        lField   = new JTextField("8, 7, 9, 8, 7, 8, 8, 6, 9, 9, 9, 7, 9, 6, 7, 7, 8," +
                " 7, 9, 7, 9, 5, 6, 6, 8, 7, 7, 8, 6, 5, 7, 7, 10, 5," +
                " 6, 5, 6, 7, 5, 6, 8");
        rField   = new JTextField("0.72363754, 0.65135398, 0.79674574, 0.69468584, 0.74442395," +
                " 0.72748567, 0.61427651, 0.65477932, 0.6619443 , 0.73156973," +
                " 0.63021967, 0.66195206, 0.69852031, 0.69098019, 0.65741823," +
                " 0.68618529, 0.75644166, 0.65312796, 0.75663052, 0.65274948," +
                " 0.72860079, 0.73703561, 0.7839319 , 0.75942205, 0.72008043," +
                " 0.71230719, 0.7315298 , 0.70147423, 0.57798889, 0.71766172," +
                " 0.71055233, 0.69250801, 0.6521953 , 0.66585477, 0.61215604," +
                " 0.67304647, 0.73438835, 0.65560526, 0.66807371, 0.65957549," +
                " 0.7917005");
        rhoField = new JTextField("0.63, 0.55, 0.56, 0.63, 0.44, 0.17, 0.51, 0.83, 0.46, 0.59, 0.93," +
                " 0.79, 0.64, 0.57, 0.97, 0.68, 0.7 , 0.26, 0.23, 0.29, 0.25, 0.19," +
                " 0.19, 0.66, 0.09, 0.8 , 0.84, 0.83, 0.14, 0.46, 0.06, 0.91, 0.92," +
                " 0.56, 0.42, 0.37, 0.93, 0.87, 0.11, 0.29, 0.82");*/

        fileField = new JTextField();
        fileField.setEditable(false);

        browseButton = new JButton("Browse MPs CSV");
        sendButton   = new JButton("Send Data");

        //Here is setting elements in correct order
        add(new JLabel("W (states per component):"));
        add(wField);

        add(new JLabel("C (capacities):"));
        add(cField);

        add(new JLabel("L (lead times):"));
        add(lField);

        add(new JLabel("R (reliabilities):"));
        add(rField);

        add(new JLabel("rho (correlations):"));
        add(rhoField);

        add(fileField);
        add(browseButton);

        add(sendButton);

        browseButton.addActionListener(e -> chooseFile());
        sendButton.addActionListener(e -> sendData());
    }

    //Implementing button for choosing a file that contains Minimal Paths
    private void chooseFile() {
        JFileChooser chooser = new JFileChooser();
        int result = chooser.showOpenDialog(this);

        if (result == JFileChooser.APPROVE_OPTION) {
            File file = chooser.getSelectedFile();
            fileField.setText(file.getAbsolutePath());
        }
    }

    //Sending data to SSVGenerator Agent. This is not agents communication yet
    private void sendData() {
        try {
            int[] W = parseIntArray(wField.getText());
            double[] C = parseDoubleArray(cField.getText());
            int[] L = parseIntArray(lField.getText());
            double[] R = parseDoubleArray(rField.getText());
            double[] rho = parseDoubleArray(rhoField.getText());

            String mpFile = fileField.getText();

            if (mpFile == null || mpFile.isEmpty()) {
                throw new IllegalArgumentException("MPs file not selected");
            }

            agent.onSendData(W, C, L, R, rho, mpFile);
            dispose();

        } catch (Exception ex) {
            JOptionPane.showMessageDialog(
                    this,
                    "Invalid input:\n" + ex.getMessage(),
                    "Error",
                    JOptionPane.ERROR_MESSAGE
            );
        }
    }

    //Methods below are used for translating text message for ints and doubles
    private int[] parseIntArray(String text) {
        String[] parts = text.split(",");
        int[] arr = new int[parts.length];
        for (int i = 0; i < parts.length; i++) {
            arr[i] = Integer.parseInt(parts[i].trim());
        }
        return arr;
    }

    private double[] parseDoubleArray(String text) {
        String[] parts = text.split(",");
        double[] arr = new double[parts.length];
        for (int i = 0; i < parts.length; i++) {
            arr[i] = Double.parseDouble(parts[i].trim());
        }
        return arr;
    }

}
