short int
linkingnumber(void) {
    int v0, v1, v2, v3, edge;
    short int l = 0;
    v0 = Vedges2[0][0];
    v1 = Vedges2[0][1];
    v2 = Vedges2[0][2];
    v3 = Vedges2[0][3];

    // First check ordertemplate vs Vedges2

    //side = 0

    edge = ordertemplate[0][1][1];

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
    edge = ordertemplate[0][1][2];
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
    edge = ordertemplate[0][1][3];
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

    //side = 1

    edge = ordertemplate[1][0][0];

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

    edge = ordertemplate[1][0][1];
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

    edge = ordertemplate[1][0][2];
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

    //Now check ordertemplate2 vs Vedges

    //side = 0

    edge = ordertemplate2[0][1][1];
    v0 = Vedges[0][0];
    v1 = Vedges[0][1];
    v2 = Vedges[0][2];
    v3 = Vedges[0][3];

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
    edge = ordertemplate2[0][1][2];
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
    edge = ordertemplate2[0][1][3];
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

    //side = 1

    edge = ordertemplate2[1][0][0];

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

    edge = ordertemplate2[1][0][1];
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

    edge = ordertemplate2[1][0][2];
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
