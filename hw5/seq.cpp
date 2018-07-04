    // interpolation_search
    auto start = chrono::high_resolution_clock::now();
    for (int i = n - 1; i >= 0; --i) {
        for (int j = 0; j <= i; ++j) {
            Node node_ij = tree[i][j];
            for (int k = 0; k <= M; ++k) {
                double A_ijk = tree[i][j].A_vec[k];
                double A_log_ijk = tree[i][j].A_log_vec[k];
                double A_u = ( (i + passing_period + 1) * A_ijk + S_t * pow(u, i + 1 - j) * pow(d, j)  ) / (i + passing_period + 2);
                double A_u_log = ( (i + passing_period + 1) * A_log_ijk + S_t * pow(u, i + 1 - j) * pow(d, j)  ) / (i + passing_period + 2);

                
                if ( abs(tree[i+1][j].A_vec[0] - tree[i+1][j].A_vec[M]) < 0.0001 ) {
                    w_u = 1;
                } else {
                    w_u = ( tree[i+1][j].A_vec[0] - A_u ) / ( tree[i+1][j].A_vec[0] - tree[i+1][j].A_vec[M]);
                }

                //w_u = ( tree[i+1][j].A_vec[0] - A_u ) / ( tree[i+1][j].A_vec[0] - tree[i+1][j].A_vec[M]);
                C_u = w_u * tree[i+1][j].C_vec[M] + (1 - w_u) * tree[i+1][j].C_vec[0];

                w_u_log = ( tree[i+1][j].A_log_vec[0] - A_u_log ) / ( tree[i+1][j].A_log_vec[0] - tree[i+1][j].A_log_vec[M]);
                C_u_log = w_u_log * tree[i+1][j].C_log_vec[M] + (1 - w_u_log) * tree[i+1][j].C_log_vec[0];

                // find A_d and k_d
                // then compute w_d, C_d
                double A_d = ( (i + passing_period + 1) * A_ijk + S_t * pow(u, i - j) * pow(d, j + 1) ) / (double) (i + passing_period + 2);
                double A_d_log = ( (i + passing_period + 1) * A_log_ijk + S_t * pow(u, i - j) * pow(d, j + 1) ) / (double) (i + passing_period + 2);
                //double A_d = ( (i + passing_period ) * A_ijk + S_t * pow(u, i - j) * pow(d, j + 1) ) / (double) (i + passing_period + 1);
                //double A_d_log = ( (i + passing_period ) * A_log_ijk + S_t * pow(u, i - j) * pow(d, j + 1) ) / (double) (i + passing_period + 1);

                if ( abs( tree[i+1][j+1].A_vec[0] - tree[i+1][j+1].A_vec[M]) < 0.0001 ) {
                    w_d = 1;
                } else {
                    w_d = ( tree[i+1][j+1].A_vec[0] - A_d ) / ( tree[i+1][j+1].A_vec[0] - tree[i+1][j+1].A_vec[M] );
                }
                //w_d = ( tree[i+1][j+1].A_vec[0] - A_d ) / ( tree[i+1][j+1].A_vec[0] - tree[i+1][j+1].A_vec[M] );

                C_d = w_d * tree[i+1][j+1].C_vec[M] + (1 - w_d) * tree[i+1][j+1].C_vec[0];
                
                w_d_log = ( tree[i+1][j+1].A_log_vec[0] - A_d_log ) / ( tree[i+1][j+1].A_log_vec[0] - tree[i+1][j+1].A_log_vec[M]);
                C_d_log = w_d_log * tree[i+1][j+1].C_log_vec[M] + (1 - w_d_log) * tree[i+1][j+1].C_log_vec[0];

                // update C(i, j, k)
                tree[i][j].C_vec[k] = (p * C_u + (1 - p) * C_d) * exp( -r * dT );
                tree[i][j].C_log_vec[k] = (p * C_u_log + (1 - p) * C_d_log) * exp( -r * dT );
                // American: max( A(i, j, k) - K, (P * C_u + (1 - P) * C_d) * e^-r * dT
                tree[i][j].C_am_vec[k] = max( A_ijk - K, tree[i][j].C_vec[k] );
                tree[i][j].C_am_log_vec[k] = max( A_log_ijk - K, tree[i][j].C_log_vec[k] );

            }
        }
    }
