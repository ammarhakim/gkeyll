#include <stdio.h>
#include <stdlib.h>

void nodal_to_modal(const double *fnodal, double *fmodal) {
    printf("Inside: %.7e %.7e %.7e %.7e\n", fnodal[0], fnodal[1], fnodal[2], fnodal[3]);
}

int main() {
    int num_basis = 4;
    int i;
    for(int j=0; j<2; j++) {
        if(1) {
            double fn_up_all[num_basis];
            for (i=0; i<num_basis; i++) fn_up_all[i] = (double)i + 1.5;
            printf("Outside: %.7e %.7e %.7e %.7e\n", fn_up_all[0], fn_up_all[1], fn_up_all[2], fn_up_all[3]);
            nodal_to_modal(fn_up_all, NULL);
        }
    }
    return 0;
}
