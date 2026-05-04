void addtowalk(int printablewalk[], int nextwalk[]) {
/* This function nextwalk to printablewalk
*/
	int i=0;
	int j=0;
	int temp1, temp2;
	int index;
	
/*
	printf("addtowalk() is called\n");
	printf("printablewalk is: ");
	while(printablewalk[i]!=0){
		printf("%d, ", printablewalk[i]);
		i++;
	}
	printf("\nnextwalk is: ");

	while(nextwalk[j]!=0){
		printf("%d, ", nextwalk[j]);
		j++;
	}
	printf("\n");
*/
	i=0;
	j=0;
	while(printablewalk[i]!=0){
		while(printablewalk[i]!=-222){
			i++;
		}
		//insert nextwalk here (printablewalk[i] is now -222)
		//shift printable walk, then add nextwalk.

		while(nextwalk[j]!=-333){	//don't add nextwalk[j] if any of the following 3 conditions are true, since it would be double counting steps.
			if(j==0){	//first entry in nextwalk (know we just entered the 2-span from the left
				j++;
			}
			else if((nextwalk[j]==5 && nextwalk[j-1]<0)){	//just entered the 2-span from the left
				j++;
			}
			else if(nextwalk[j]==6 && nextwalk[j+1]==-333){	//just left the 2-span to the left
				j++;
			}
			else{	//if the above 3 conditions fail, then we need to add nextwalk[j] to printablewalk
				index=i;
//				printf("RUNNING. i=%d, index=%d, j=%d\n", i, index, j);
				//shift printablewalk one over
				temp1=printablewalk[index];
				while(printablewalk[index]!=0){
					temp2=printablewalk[index+1];
					printablewalk[index+1]=temp1;
					temp1=temp2;
					index++;
				}
				//add nextwalk component
				printablewalk[i]=nextwalk[j];
				j++;
				i++;
			}
		}
//		printf("after finished adding:, i=%d, j=%d, index=%d, \n", i, j, index);
		//now delete all negative numbers starting from printablewalk[i] onwards
		while(printablewalk[i]<0){
//			printf("deleting printablewalk[%d]=%d\n", i, printablewalk[i]);
			index=i;
			while(printablewalk[index]!=0){
				printablewalk[index] = printablewalk[index+1];
				index++;
			}
		}

		//advance j until nextwalk[j] >=0
		while(nextwalk[j]<0){
			j++;
		}
//		printf("advanced j. nextwalk[%d]=%d\n", j, nextwalk[j]);

		if(nextwalk[j]==0){	//then no more needed to add
			break;
		}

	}

/*
	printf("AFTER INSERTION: ");
	i=0;
	j=0;
	printf("printablewalk is: ");
	while(printablewalk[i]!=0){
		printf("%d, ", printablewalk[i]);
		i++;
	}
	printf("\n");
*/
/*
	printf("\nnextwalk is: ");

	while(nextwalk[j]!=0){
		printf("%d, ", nextwalk[j]);
		j++;
	}
	printf("\n");
*/





}
