import edu.princeton.cs.algs4.Picture;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.Stopwatch;

import java.awt.*;

/**
 * Created by xgy on 21/06/16.
 */
public class SeamCarver {
    private Picture pic;
    private Color[][] picArr;
    private boolean transposed;
    private int trueWidth;
    private int trueHeight;

    /**
     * Constructor for SeamCarver class.
     * It reads a Picture object and save it to a 2D array.
     *
     * @param picture
     */
    public SeamCarver(Picture picture) {
        this.pic = picture;
        this.picArr = new Color[picture.height()][picture.width()];

        for (int row = 0; row < picArr.length; row++) {
            for (int col = 0; col < picArr[0].length; col++) {
                picArr[row][col] = picture.get(col, row);
            }
        }
        transposed = false;
        trueWidth = picArr[0].length;
        trueHeight = picArr.length;
    }

    /**
     * Get current picture to client.
     *
     * @return Picture type
     */
    public Picture picture() {
        this.pic = arrToPic();
        return this.pic;
    }

    /**
     * Transform from 2D array to Picture object.
     *
     * @return Picture object
     */
    private Picture arrToPic() {
        if (transposed) deTranspose();
        Picture p = new Picture(trueWidth, trueHeight);
        for (int row = 0; row < p.height(); row++) {
            for (int col = 0; col < p.width(); col++) {
                p.set(col, row, picArr[row][col]);
            }
        }
        return p;
    }

    /**
     * Return current width, not array size
     *
     * @return int value of width
     */
    public int width() {
        if (transposed) return trueHeight;
        else return trueWidth;
    }

    /**
     * return current height, not array size
     *
     * @return int
     */
    public int height() {
        if (transposed) return trueWidth;
        else return trueHeight;
    }

    /**
     * public api to get cell's energy, it translate coordinates to current coordinates.
     *
     * @param col
     * @param row
     * @return double: energy
     */
    public double energy(int col, int row) {
        if (transposed) {
            return getEnergy(trueWidth - 1 - row, col);
        } else {
            return getEnergy(col, row);
        }
    }

    /**
     * Private helper function to get energy.
     *
     * @param col
     * @param row
     * @return
     */
    private double getEnergy(int col, int row) {
        if (col == 0 || col == trueWidth - 1 || row == 0 || row == trueHeight - 1) {
            return 1000;
        }
        double deltaX = deltaX2(col, row);
        double deltaY = deltaY2(col, row);
        return Math.sqrt(deltaX + deltaY);
    }


    /**
     * Calculate horizontal energy from left and right cells.
     *
     * @param col
     * @param row
     * @return
     */
    private double deltaX2(int col, int row) {
        double redDiff = picArr[row][col - 1].getRed() - picArr[row][col + 1].getRed();
        double greenDiff = picArr[row][col - 1].getGreen() - picArr[row][col + 1].getGreen();
        double blueDiff = picArr[row][col - 1].getBlue() - picArr[row][col + 1].getBlue();

        double res = Math.pow(redDiff, 2) + Math.pow(greenDiff, 2) + Math.pow(blueDiff, 2);
        if (Double.isNaN(res)) return 0.0;
        else return res;
    }

    /**
     * Calculate vertical energy from up and down cells.
     *
     * @param col
     * @param row
     * @return
     */
    private double deltaY2(int col, int row) {
        double redDiff = picArr[row - 1][col].getRed() - picArr[row + 1][col].getRed();
        double greenDiff = picArr[row - 1][col].getGreen() - picArr[row + 1][col].getGreen();
        double blueDiff = picArr[row - 1][col].getBlue() - picArr[row + 1][col].getBlue();

        double res = Math.pow(redDiff, 2) + Math.pow(greenDiff, 2) + Math.pow(blueDiff, 2);
        if (Double.isNaN(res)) return 0.0;
        else return res;
    }

    /**
     * @param col
     * @param row
     * @param energy  2D array recording energy of each cell
     * @param distTo  previous cell pointing to this cell, which make a shortest path(sum of energy) to this cell.
     * @param pathSum path sum (sum of energy)
     */
    private void verticalRelax(int col, int row, double[][] energy, int[][][] distTo, double[][] pathSum) {
        /* look into three cell on the next row */
        double curSumEnergy = pathSum[col][row];

        /* down-left cell */
        if (col - 1 >= 0) {
            double targetSumEnergy = pathSum[col - 1][row + 1];
            if (targetSumEnergy == 0.0 || targetSumEnergy > energy[col - 1][row + 1] + curSumEnergy) {
                pathSum[col - 1][row + 1] = energy[col - 1][row + 1] + curSumEnergy;
                distTo[col - 1][row + 1][0] = col;
                distTo[col - 1][row + 1][1] = row;
            }
        }
        /* down cell */
        double targetSumEnergyMid = pathSum[col][row + 1];
        if (targetSumEnergyMid == 0.0 || targetSumEnergyMid > energy[col][row + 1] + curSumEnergy) {
            pathSum[col][row + 1] = energy[col][row + 1] + curSumEnergy;
            distTo[col][row + 1][0] = col;
            distTo[col][row + 1][1] = row;
        }

        /* down-right cell */
        if (col + 1 < trueWidth) {
            double targetSumEnergy = pathSum[col + 1][row + 1];
            if (targetSumEnergy == 0.0 || targetSumEnergy > energy[col + 1][row + 1] + curSumEnergy) {
                pathSum[col + 1][row + 1] = energy[col + 1][row + 1] + curSumEnergy;
                distTo[col + 1][row + 1][0] = col;
                distTo[col + 1][row + 1][1] = row;
            }
        }
    }


    /**
     * Find shortest path from up to bottom
     *
     * @return an array of index of col
     */
    private int[] findSeam() {
        int[] res = new int[trueHeight];
        double[][] energy = new double[trueWidth][trueHeight];
        double[][] pathSum = new double[trueWidth][trueHeight];

        /**
         * (cur.x, cur.y) -> [from.x, from.y]
         */
        int[][][] distTo = new int[trueWidth][trueHeight][2];

        /* calculate energy */
        for (int row = 0; row < trueHeight; row++) {
            for (int col = 0; col < trueWidth; col++) {
                energy[col][row] = getEnergy(col, row);
            }
        }

        /* relax edges updating shortest path */
        for (int row = 0; row < trueHeight - 1; row++) {
            for (int col = 0; col < trueWidth; col++) {
                verticalRelax(col, row, energy, distTo, pathSum);
            }
        }


        double minPathSum = Double.MAX_VALUE;
        int minIndex = 0;
        for (int col = 0; col < trueWidth; col++) {
            if (pathSum[col][trueHeight - 1] < minPathSum) {
                minPathSum = pathSum[col][trueHeight - 1];
                minIndex = col;
            }
        }

        res[res.length - 1] = minIndex;
        int count = res.length - 1;
        while (count > 0) {
            /* trace back from bottom to up */
            int fromCol = distTo[minIndex][count][0];
            res[count - 1] = fromCol;
            minIndex = fromCol;
            --count;
        }

        return res;
    }

    /**
     * Transpose current pic array.
     */
    private void transpose() {
        Color[][] transposedArr = new Color[trueWidth][trueHeight];
        for (int col = 0; col < trueWidth; col++) {
            for (int row = 0; row < trueHeight; row++) {
                transposedArr[col][trueHeight - 1 - row] = picArr[row][col];
            }
        }
        picArr = transposedArr;
        transposed = !transposed;
        int tmp = trueHeight;
        trueHeight = trueWidth;
        trueWidth = tmp;
    }

    /**
     * Back to original angle
     */
    private void deTranspose() {
        Color[][] deTransposedArr = new Color[trueWidth][trueHeight];
        for (int col = 0; col < trueWidth; col++) {
            for (int row = 0; row < trueHeight; row++) {
                deTransposedArr[trueWidth - 1 - col][row] = picArr[row][col];
            }
        }
        picArr = deTransposedArr;
        transposed = !transposed;
        int tmp = trueHeight;
        trueHeight = trueWidth;
        trueWidth = tmp;
    }

    /**
     * @return
     */
    public int[] findVerticalSeam() {
        if (transposed) {
            deTranspose();
        }
        return findSeam();
    }


    /**
     * Transpose it and use the same way as finding vertical seam
     *
     * @return
     */
    public int[] findHorizontalSeam() {
        if (!transposed) {
            transpose();
        }
        int[] res = findSeam();
        for (int i = 0; i < res.length; i++) {
            res[i] = trueWidth - 1 - res[i];
        }
        return res;
    }

    public void removeHorizontalSeam(int[] seam) {
        if (!transposed) {
            transpose();
        }
        for (int i = 0; i < seam.length; i++) {
            seam[i] = trueWidth - 1 - seam[i];
        }

        for (int row = 0; row < seam.length; row++) {
            System.arraycopy(picArr[row], seam[row] + 1, picArr[row], seam[row], trueWidth - seam[row] - 1);
        }
        trueWidth--;

    }

    /**
     * Given seam calculated previously, cut the seam
     *
     * @param seam
     */
    public void removeVerticalSeam(int[] seam) {
        if (transposed) {
            deTranspose();
        }

        for (int row = 0; row < seam.length; row++) {
            System.arraycopy(picArr[row], seam[row] + 1, picArr[row], seam[row], trueWidth - seam[row] - 1);
        }
        trueWidth--;
    }


    public static void main(String[] args) {
        if (args.length != 3) {
            StdOut.println("Usage:\njava ResizeDemo [image filename] [num cols to remove] [num rows to remove]");
            return;
        }

        Picture inputImg = new Picture(args[0]);
        int removeColumns = Integer.parseInt(args[1]);
        int removeRows = Integer.parseInt(args[2]);

        StdOut.printf("image is %d columns by %d rows\n", inputImg.width(), inputImg.height());
        SeamCarver sc = new SeamCarver(inputImg);

        Stopwatch sw = new Stopwatch();

        for (int i = 0; i < removeRows; i++) {
            int[] horizontalSeam = sc.findHorizontalSeam();
            sc.removeHorizontalSeam(horizontalSeam);
        }

        for (int i = 0; i < removeColumns; i++) {
            int[] verticalSeam = sc.findVerticalSeam();
            sc.removeVerticalSeam(verticalSeam);
        }
        Picture outputImg = sc.picture();

        StdOut.printf("new image size is %d columns by %d rows\n", sc.width(), sc.height());

        StdOut.println("Resizing time: " + sw.elapsedTime() + " seconds.");
        inputImg.show();
        outputImg.show();
//        Picture picture = new Picture("/home/xgy/IdeaProjects/SeamCarver/src/6x5.png");
//        StdOut.printf("image is %d pixels wide by %d pixels high.\n", picture.width(), picture.height());
//
//        SeamCarver sc = new SeamCarver(picture);
//
//        StdOut.printf("Printing energy calculated for each pixel.\n");
//
//        for (int row = 0; row < sc.height(); row++) {
//            for (int col = 0; col < sc.width(); col++)
//                StdOut.printf("%9.0f ", sc.energy(col, row));
//            StdOut.println();
//        }
//        int[] res = sc.findVerticalSeam();
//        System.out.println("");
//        for (int each : res) {
//            System.out.println(each);
//        }
//        int[] res2 = sc.findHorizontalSeam();
//        System.out.println("");
//        for (int each : res2) {
//            System.out.println(each);
//        }
//
//        System.out.println("\n" + sc.height());
////        sc.removeHorizontalSeam();
//        sc.removeVerticalSeam(sc.findVerticalSeam());
//        System.out.println(sc.height());
////        assert (picture.height() == sc.picture().height());
//        for (int row = 0; row < sc.height(); row++) {
//            for (int col = 0; col < sc.width(); col++)
//                StdOut.printf("%9.0f ", sc.energy(col, row));
//            StdOut.println();
//        }

    }
}
