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
        setTitle("SSV Generator");
        setSize(600, 400);
        setLocationRelativeTo(null);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        setLayout(new GridLayout(7, 2, 5, 5));

        // Default values
        wField   = new JTextField("4,3,2,3,2");
        cField   = new JTextField("10,15,25,15,20");
        lField   = new JTextField("5,7,6,5,8");
        rField   = new JTextField("0.7,0.65,0.67,0.71,0.75");
        rhoField = new JTextField("0.1,0.3,0.5,0.7,0.9");

        fileField = new JTextField();
        fileField.setEditable(false);

        browseButton = new JButton("Browse MPs CSV");
        sendButton   = new JButton("Send Data");

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

    private void chooseFile() {
        JFileChooser chooser = new JFileChooser();
        int result = chooser.showOpenDialog(this);

        if (result == JFileChooser.APPROVE_OPTION) {
            File file = chooser.getSelectedFile();
            fileField.setText(file.getAbsolutePath());
        }
    }

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
