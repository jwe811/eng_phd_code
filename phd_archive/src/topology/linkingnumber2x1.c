short int
linkingnumber(void)
{
    int a,b;
    short int l=0;

    a = ordertemplate[0][1][1];
    b = Vedges2[0][0];

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
    a = ordertemplate[0][1][2];
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
    b = Vedges2[0][1];
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
    a = ordertemplate[1][0][0];
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
    b = Vedges2[0][2];
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
    a = ordertemplate[1][0][1];
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


    a = ordertemplate2[0][1][1];
    b = Vedges[0][0];

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
    a = ordertemplate2[0][1][2];
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
    b = Vedges[0][1];
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
    a = ordertemplate2[1][0][0];
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
    b = Vedges[0][2];
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
    a = ordertemplate2[1][0][1];
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
