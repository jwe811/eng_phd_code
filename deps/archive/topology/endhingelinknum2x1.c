
short int
endlinkingnumberright(void){
    int a,b;
    short int l=0;

    a = endtemplate1[1][1];
    b = Vedgesend2[0][0];

    if (a!=0 && b!=0){
        if (a%2 == 0){
            if (b<0){
                l+=1;
            }
            else{
                l-=1;
            }
        }
        else{
            if(b<0){
                l-=1;
            }
            else{
                l+=1;
            }
        }
    }
    a = endtemplate1[1][2];
    if (a!=0 && b!=0){
        if (a%2 == 0){
            if (b<0){
                l+=1;
            }
            else{
                l-=1;
            }
        }
        else{
            if(b<0){
                l-=1;
            }
            else{
                l+=1;
            }
        }
    }
    b = Vedgesend2[0][1];
    if (a!=0 && b!=0){
        if (a%2 == 0){
            if (b<0){
                l+=1;
            }
            else{
                l-=1;
            }
        }
        else {
            if (b < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        }
    }

    a = endtemplate2[1][1];
    b = Vedgesend[0][0];
    if (a!=0 && b!=0){
        if (a%2 == 0){
            if (b<0){
                l+=1;
            }
            else{
                l-=1;
            }
        }
        else{
            if(b<0){
                l-=1;
            }
            else{
                l+=1;
            }
        }
    }
    a = endtemplate2[1][2];
    if (a!=0 && b!=0){
        if (a%2 == 0){
            if (b<0){
                l+=1;
            }
            else{
                l-=1;
            }
        }
        else{
            if(b<0){
                l-=1;
            }
            else{
                l+=1;
            }
        }
    }
    b = Vedgesend[0][1];
    if (a!=0 && b!=0){
        if (a%2 == 0){
            if (b<0){
                l+=1;
            }
            else{
                l-=1;
            }
        }
        else {
            if (b < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        }
    }
    return l;

}
short int
endlinkingnumberleft(void){
    int a, b;
    short int l=0;

    a = endtemplate3[0][0];
    b = Vedgesend2[0][1];

    if (a!=0 && b!=0){
        if (a%2 == 0){
            if (b<0){
                l-=1;
            }
            else{
                l+=1;
            }
        }
        else{
            if(b<0){
                l+=1;
            }
            else{
                l-=1;
            }
        }
    }

    b = Vedgesend2[0][2];
    if (a!=0 && b!=0){
        if (a%2 == 0){
            if (b<0){
                l-=1;
            }
            else {

                l += 1;
            }
        }
        else{

            if(b<0){
                l+=1;
            }
            else{
                l-=1;
            }
        }
    }

    a = endtemplate3[0][1];
    if (a!=0 && b!=0) {
        if (a % 2 == 0) {
            if (b < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        } else {
            if (b < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        }
    }

    a = endtemplate4[0][0];
    b = Vedgesend[0][1];

    if (a!=0 && b!=0){
        if (a%2 == 0){
            if (b<0){
                l-=1;
            }
            else{
                l+=1;
            }
        }
        else{
            if(b<0){
                l+=1;
            }
            else{
                l-=1;
            }
        }
    }

    b = Vedgesend[0][2];
    if (a!=0 && b!=0){
        if (a%2 == 0){
            if (b<0){
                l-=1;
            }
            else {

                l += 1;
            }
        }
        else{

            if(b<0){
                l+=1;
            }
            else{
                l-=1;
            }
        }
    }

    a = endtemplate4[0][1];
    if (a!=0 && b!=0) {
        if (a % 2 == 0) {
            if (b < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        } else {
            if (b < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        }
    }
    return l;
}