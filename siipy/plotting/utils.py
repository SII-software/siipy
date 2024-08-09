def cen(f, zoom=1):
    n = f.shape[0]
    m = n // (zoom * 2)
    return f[n // 2 - m:n // 2 + m, n // 2 - m:n // 2 + m]