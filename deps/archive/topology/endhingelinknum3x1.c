
short int
endlinkingnumberright(void){
    int v0,v1,v2,v3, edge;
    short int l=0;

    edge = endtemplate1[1][1];
    v0 = Vedgesend2[0][0];
    v1 = Vedgesend2[0][1];
    v2 = Vedgesend2[0][2];
    v3 = Vedgesend2[0][3];

    //First check endtemplate1 vs Vedgesend2

    edge = endtemplate1[1][1];

    if (edge != 0 && v0 != 0) {
        if (edge % 2 == 0) {
            if (v0 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        } else {
            if (v0 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        }
    }
    edge = endtemplate1[1][2];
    if (edge != 0 && v0 != 0) {
        if (edge % 2 == 0) {
            if (v0 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        } else {
            if (v0 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        }
    }

    if (edge != 0 && v1 != 0) {
        if (edge % 2 == 0) {
            if (v1 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        } else {
            if (v1 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        }
    }
    edge = endtemplate1[1][3];
    if (edge != 0 && v0 != 0) {
        if (edge % 2 == 0) {
            if (v0 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        } else {
            if (v0 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        }
    }
    if (edge != 0 && v1 != 0) {
        if (edge % 2 == 0) {
            if (v1 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        } else {
            if (v1 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        }
    }
    if (edge != 0 && v2 != 0) {
        if (edge % 2 == 0) {
            if (v2 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        } else {
            if (v2 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        }
    }

    //Now check endtemplate2 vs Vedgesend
    v0 = Vedgesend[0][0];
    v1 = Vedgesend[0][1];
    v2 = Vedgesend[0][2];
    v3 = Vedgesend[0][3];

    edge = endtemplate2[1][1];

    if (edge != 0 && v0 != 0) {
        if (edge % 2 == 0) {
            if (v0 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        } else {
            if (v0 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        }
    }
    edge = endtemplate2[1][2];
    if (edge != 0 && v0 != 0) {
        if (edge % 2 == 0) {
            if (v0 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        } else {
            if (v0 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        }
    }

    if (edge != 0 && v1 != 0) {
        if (edge % 2 == 0) {
            if (v1 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        } else {
            if (v1 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        }
    }
    edge = endtemplate2[1][3];
    if (edge != 0 && v0 != 0) {
        if (edge % 2 == 0) {
            if (v0 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        } else {
            if (v0 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        }
    }
    if (edge != 0 && v1 != 0) {
        if (edge % 2 == 0) {
            if (v1 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        } else {
            if (v1 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        }
    }
    if (edge != 0 && v2 != 0) {
        if (edge % 2 == 0) {
            if (v2 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        } else {
            if (v2 < 0) {
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
    int v0, v1, v2, v3, edge;
    short int l=0;

    v0 = Vedgesend2[0][0];
    v1 = Vedgesend2[0][1];
    v2 = Vedgesend2[0][2];
    v3 = Vedgesend2[0][3];

    //First check endtemplate3 vs Vedgesend2

    edge = endtemplate3[0][0];

    if (edge != 0 && v1 != 0) {
        if (edge % 2 == 0) {
            if (v1 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        } else {
            if (v1 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        }
    }

    if (edge != 0 && v2 != 0) {
        if (edge % 2 == 0) {
            if (v2 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        } else {
            if (v2 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        }
    }

    if (edge != 0 && v3 != 0) {
        if (edge % 2 == 0) {
            if (v3 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        } else {
            if (v3 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        }
    }

    edge = endtemplate3[0][1];
    if (edge != 0 && v2 != 0) {
        if (edge % 2 == 0) {
            if (v2 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        } else {
            if (v2 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        }
    }

    if (edge != 0 && v3 != 0) {
        if (edge % 2 == 0) {
            if (v3 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        } else {
            if (v3 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        }
    }

    edge = endtemplate3[0][2];
    if (edge != 0 && v3 != 0) {
        if (edge % 2 == 0) {
            if (v3 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        } else {
            if (v3 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        }
    }

    //Now check endtemplate4 vs Vedgesend
    v0 = Vedgesend[0][0];
    v1 = Vedgesend[0][1];
    v2 = Vedgesend[0][2];
    v3 = Vedgesend[0][3];

    edge = endtemplate4[0][0];

    if (edge != 0 && v1 != 0) {
        if (edge % 2 == 0) {
            if (v1 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        } else {
            if (v1 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        }
    }

    if (edge != 0 && v2 != 0) {
        if (edge % 2 == 0) {
            if (v2 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        } else {
            if (v2 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        }
    }

    if (edge != 0 && v3 != 0) {
        if (edge % 2 == 0) {
            if (v3 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        } else {
            if (v3 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        }
    }

    edge = endtemplate4[0][1];
    if (edge != 0 && v2 != 0) {
        if (edge % 2 == 0) {
            if (v2 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        } else {
            if (v2 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        }
    }

    if (edge != 0 && v3 != 0) {
        if (edge % 2 == 0) {
            if (v3 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        } else {
            if (v3 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        }
    }

    edge = endtemplate4[0][2];
    if (edge != 0 && v3 != 0) {
        if (edge % 2 == 0) {
            if (v3 < 0) {
                l -= 1;
            } else {
                l += 1;
            }
        } else {
            if (v3 < 0) {
                l += 1;
            } else {
                l -= 1;
            }
        }
    }


    return l;
}