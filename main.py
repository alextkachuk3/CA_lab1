import numpy as np


def poly_to_str(poly):
    str_poly = ''
    if len(poly) != 0:
        for i in range(len(poly) - 1, 0, -1):
            if poly[len(poly) - i - 1] >= 0:
                if i == len(poly) - 1:
                    if poly[len(poly) - i - 1] == 1:
                        str_poly = str_poly + 'x'
                    else:
                        str_poly = str_poly + '{}x'.format(poly[len(poly) - i - 1])
                    if i != 1:
                        str_poly = str_poly + '^{}'.format(i)
                else:
                    if poly[len(poly) - i - 1] == 1:
                        str_poly = str_poly + '+x'
                    else:
                        str_poly = str_poly + '+{}x'.format(poly[len(poly) - i - 1])
                    if i != 1:
                        str_poly = str_poly + '^{}'.format(i)
            else:
                if poly[len(poly) - i - 1] == -1:
                    str_poly = str_poly + '-x'
                else:
                    str_poly = str_poly + '{}x'.format(poly[len(poly) - i - 1])
                if i != 1:
                    str_poly = str_poly + '^{}'.format(i)
        if poly[len(poly) - 1] >= 0:
            str_poly = str_poly + '+{}'.format(poly[-1])
        else:
            str_poly = str_poly + '{}'.format(poly[-1])
    return str_poly


def poly_list_to_str(poly_list):
    res = '('
    for i in range(0, len(poly_list)):
        res = res + poly_to_str(poly_list[i]) + ')'
        if i != len(poly_list) - 1:
            res = res + '('
    return res


def kroneker_factorization(poly):
    res = []
    poly_f = kroneker(poly)
    if poly_f:
        poly_div = np.polydiv(poly, poly_f)[0].astype(int).tolist()
        res.append(poly_f)
    while poly_f is not False:
        poly_f = kroneker(poly_div)
        if poly_f:
            poly_div = np.polydiv(poly_div, poly_f)[0].astype(int).tolist()
            res.append(poly_f)
        else:
            res.append(poly_div)
    return res


def kroneker(poly):
    print()
    if len(poly) <= 2:
        print('Polynom ({}) is irreducible'.format(poly_to_str(poly)))
        return False

    print('Factorization {}'.format(poly_to_str(poly)))
    m = np.fix((len(poly) - 1) / 2).astype(int)
    p_i = []
    x = []
    if m > 0:
        for i in range(0, m + 1):
            x.append(i)
            cur_p_i = calc_poly(poly, i)
            p_i.append(cur_p_i)
            output = "i={} P(i)={}".format(i, cur_p_i)
            if cur_p_i == 0:
                poly_div = [1, -i]
                output = output + ' => ({}) is divisor'.format(poly_to_str(poly_div))
                print(output)
                return poly_div
            print(output)

    p_i_div = []
    for cur_p_i in p_i:
        p_i_div.append(divisors(cur_p_i))

    print()
    print('Divisors{}'.format(p_i_div))
    p_i_div_cartesian = cartesian(p_i_div)
    print('M={}'.format(p_i_div_cartesian))
    print('M size {}'.format(len(p_i_div_cartesian)))
    print()

    for i in range(0, len(p_i_div_cartesian)):
        output = '#{} X={} Y={} '.format(i + 1, x, p_i_div_cartesian[i])
        inter_poly = interpol(x, p_i_div_cartesian[i])
        output = output + 'g(x)={} '.format(poly_to_str(inter_poly))
        print(output)
        if len(inter_poly) > 1:
            poly_div = np.polydiv(poly, inter_poly)
            output = 'P(x) / g(x) = {}, remainder = {}'.format(poly_to_str(poly_div[0].astype(float).tolist()),
                                                               poly_div[1][0])
            if poly_div[1][0] == 0.0:
                output = output + ' => ({}) is divisor'.format(poly_to_str(inter_poly))
                print(output)
                return inter_poly

            print(output)
        else:
            print('Polynom power is not positive')

    print('Polynom ({}) is irreducible'.format(poly_to_str(poly)))
    return False


def calc_poly(poly, x):
    res = 0
    for i in range(len(poly) - 1, -1, -1):
        res = res + ((x ** i) * poly[len(poly) - i - 1])
    return res


def cartesian_product(x, y):
    result = []
    for i in range(0, len(x)):
        for j in range(0, len(y)):
            if type(x[i]) != list:
                x[i] = [x[i]]
            temp = [num for num in x[i]]
            temp.append(y[j])
            result.append(temp)
    return result


def cartesian(sets_list):
    n = len(sets_list)
    result = sets_list[0]
    for i in range(1, n):
        result = cartesian_product(result, sets_list[i])
    return result


def divisors(x):
    res = []
    i = 1
    if x < 0:
        x = -x
    while i <= x:
        if x % i == 0:
            res.append(i)
            res.append(-i)
        i = i + 1
    return res


def interpol(x, y):
    numerator = [0]
    denominator = 1
    for i in range(0, len(y)):
        cur_numerator = [y[i]]
        cur_denominator = 1
        for j in range(0, len(x)):
            if i != j:
                cur_numerator = np.polymul(cur_numerator, [1, -x[j]])
                cur_denominator = cur_denominator * (x[i] - x[j])
        numerator = np.polyadd(np.polymul(numerator, cur_denominator), np.polymul(cur_numerator, denominator))
        denominator = denominator * cur_denominator
    res = np.polydiv(numerator, denominator)[0].astype(int).tolist()

    i = 0
    while i < len(res):
        if res[i] == 0:
            res.pop(i)
        else:
            break
    return res


if __name__ == '__main__':
    poly1 = [6, 23, -9, -92, -60]
    fact1 = kroneker_factorization(poly1)
    print()
    print('{} = {}'.format(poly_to_str(poly1), poly_list_to_str(fact1)))

    poly2 = [2, -25, 93, -90]
    fact2 = kroneker_factorization(poly2)
    print()
    print('{} = {}'.format(poly_to_str(poly2), poly_list_to_str(fact2)))
