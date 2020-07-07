def tzf_wheel(n):
    R.<s> = QQ['s']
    
    # cycle through some base cases ...
    if (n == 0):
        return 1
    elif (n == 1):
        return 1/(1+s)

    # prepare the fans for recursion...
    compute_fans(n-1)

    single_sum = 0
    for k in range(1, n):
        single_sum += ((-s)/(s+1))^(n-k-1) * fan_tzfs[k]
        
    non_fans = 0
    for r in range(0, n-1):
        non_fans += binomial(n, r) * (n-1-r) * (-1)^(n-1-r) * 1/(s+1)^r

    # Brute force add cyclic matroid final case
    non_fans += ((-s)^n-1-n*(-s-1)) / ((n*s+n-1) * (s+1)^n)

    total = n * single_sum + non_fans

    return 1/(2*n*s + n) * total
