int binary_search ( vector < double > &A, double &key, int left, int right ) {
    int mid = left + ( right - left ) / 2 ;
    if ( abs ( A[mid] - key ) < 0.0000001 ) return mid ;
    else if ( A[mid] > key ) return binary_search ( A, key, mid + 1, right ) ;
    else if ( A[mid] < key ) {
        if ( abs ( A[mid-1] - key ) < 0.0000001 ) return mid-1 ;
        else if ( A[mid-1] > key ) return mid ;
        else return binary_search ( A, key, left, mid - 1 ) ;
    }
}

int linear_interpolation ( vector < double > &A, double &key, int left, int right ) {
    int mid ;
    if ( left == right ) mid = left ;
    else if ( ( A[left] - A[right] ) < 0.0000001 ) return left ;
    else mid = left + ( right - left ) * ( key - A[left] ) / ( A[right] - A[left] ) ;
    //mid = left + ( right - left ) * ( key - A[left] ) / ( A[right] - A[left] ) ;
    if ( abs ( A[mid] - key ) < 0.0000001 ) return mid ;
    else if ( A[mid] > key ) return linear_interpolation ( A, key, mid + 1, mid + 1 ) ;
    else if ( A[mid] < key ) {
        if ( abs ( A[mid-1] - key ) < 0.0000001 ) return mid-1 ;
        else if ( A[mid-1] > key ) return mid ;
        else return linear_interpolation ( A, key, mid - 1, mid - 1 ) ;
    }
}
