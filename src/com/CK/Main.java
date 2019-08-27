package com.CK;

import java.util.Arrays;

public class Main {

    public static void main(String[] args) {
        String[] A = {"catg", "ctaagt", "gcta", "ttca", "atgcatc"};
        new Solution2().shortestSuperstring(A);
    }
}

//TSP
class Solution {
    int n = 4;
    int[][] dist = {
            {0, 20, 42, 25},
            {20, 0, 30, 34},
            {42, 30, 0, 10},
            {25, 34, 10, 0}
    };
    int VISITED_ALL = (1 << n) - 1;
    int[][] dp = new int[(1 << n)][n];

    int tsp(int mask, int pos) {
        if (mask == VISITED_ALL) {
            return dist[pos][0];
        }
        if (dp[mask][pos] != -1) {
            return dp[mask][pos];
        }
        //Now from current node, we will try to go to every other node and take the min ans
        int ans = Integer.MAX_VALUE;
        //Visit all the unvisited cities and take the best route
        for (int city = 0; city < n; city++) {
            if ((mask & (1 << city)) == 0) {
                int newAns = dist[pos][city] + tsp(mask | (1 << city), city);
                ans = Math.min(ans, newAns);
            }
        }
        return dp[mask][pos] = ans;
    }

    void main() {
        /* init the dp array */
        for (int i = 0; i < (1 << n); i++) {
            for (int j = 0; j < n; j++) {
                dp[i][j] = -1;
            }
        }
        System.out.println("Travelling Saleman Distance is " + tsp(1, 0));
    }
}

class Solution2 {
    public String shortestSuperstring(String[] A) {
        int N = A.length;
        // Populate overlaps
        int[][] overlaps = new int[N][N];
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                if (i != j) {
                    int m = Math.min(A[i].length(), A[j].length());
                    for (int k = m; k >= 0; --k)
                        if (A[i].endsWith(A[j].substring(0, k))) {
                            overlaps[i][j] = k;
                            break;
                        }
                }

        // dp[mask][i] = most overlap with mask, ending with ith element
        int[][] dp = new int[1 << N][N];
        int[][] parent = new int[1 << N][N];
        for (int mask = 0; mask < (1 << N); ++mask) {
            Arrays.fill(parent[mask], -1);
            for (int bit = 0; bit < N; ++bit)
                if ((mask & (1 << bit)) > 0) {
                    // Let's try to find dp[mask][bit].  Previously, we had
                    // a collection of items represented by pmask.
                    int pmask = mask ^ (1 << bit);
                    if (pmask == 0) continue;
                    for (int i = 0; i < N; ++i)
                        if ((pmask & (1 << i)) > 0) {
                            // For each bit i in pmask, calculate the value
                            // if we ended with word i, then added word 'bit'.
                            int val = dp[pmask][i] + overlaps[i][bit];
                            if (val > dp[mask][bit]) {
                                dp[mask][bit] = val;
                                parent[mask][bit] = i;
                            }
                        }
                }
        }

        // # Answer will have length sum(len(A[i]) for i) - max(dp[-1])
        // Reconstruct answer, first as a sequence 'perm' representing
        // the indices of each word from left to right.

        int[] perm = new int[N];
        boolean[] seen = new boolean[N];
        int t = 0;
        int mask = (1 << N) - 1;

        // p: the last element of perm (last word written left to right)
        int p = 0;
        for (int j = 0; j < N; ++j)
            if (dp[(1 << N) - 1][j] > dp[(1 << N) - 1][p])
                p = j;

        // Follow parents down backwards path that retains maximum overlap
        while (p != -1) {
            perm[t++] = p;
            seen[p] = true;
            int p2 = parent[mask][p];
            mask ^= 1 << p;
            p = p2;
        }

        // Reverse perm
        for (int i = 0; i < t / 2; ++i) {
            int v = perm[i];
            perm[i] = perm[t - 1 - i];
            perm[t - 1 - i] = v;
        }

        // Fill in remaining words not yet added
        for (int i = 0; i < N; ++i)
            if (!seen[i])
                perm[t++] = i;

        // Reconstruct final answer given perm
        StringBuilder ans = new StringBuilder(A[perm[0]]);
        for (int i = 1; i < N; ++i) {
            int overlap = overlaps[perm[i - 1]][perm[i]];
            ans.append(A[perm[i]].substring(overlap));
        }

        return ans.toString();
    }
}
