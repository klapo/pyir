def undistortUV(Ud, Vd, ip):
    '''
    [U V] = undistortUV(Ud,Vd,ip)
    ip - [lx ly uWidth vWidth u0 v0 d1 d2]
    '''
    # separate parameters
    lu = ip[0]
    lv = ip[1]
    Uo = ip[4]
    Vo = ip[5]
    d1 = ip[6]
    d2 = ip[7]
    if len(ip) > 8:
        d3 = ip[8]

    # calculate radius
    r = np.sqrt((lu * Ud - Uo)**2 + (lv * Vd - Vo)**2)

    # calculate scaling factor
    if len(ip) > 8:
        s = 1 + d1 * r**2 + d3 * r + d2
    else:
        s = 1 + d1 * r**2 + d2

    # undistort UV
    U = ((lu * Ud - Uo) * s + Uo) / lu
    V = ((lv * Vd - Vo) * s + Vo) / lv

    return(U, V)
