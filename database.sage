from urllib.request import urlopen

url_base = 'http://www-imai.is.s.u-tokyo.ac.jp/~ymatsu/matroid/database/'

rank_1 = [None, 'allr1n01.txt', 'allr1n02.txt', 'allr1n03.txt', 'allr1n04.txt',
    'allr1n05.txt', 'allr1n06.txt', 'allr1n07.txt', 'allr1n08.txt', 'allr1n09.txt',
    'allr1n10.txt', 'allr1n11.txt', 'allr1n12.txt']

rank_2 = [None, None, 'allr2n02.txt', 'allr2n03.txt', 'allr2n04.txt',
    'allr2n05.txt', 'allr2n06.txt', 'allr2n07.txt', 'allr2n08.txt', 'allr2n09.txt',
    'allr2n10.txt', 'allr2n11.txt', 'allr2n12.txt']

rank_3 = [None, None, None, 'allr3n03.txt', 'allr3n04.txt',
    'allr3n05.txt', 'allr3n06.txt', 'allr3n07.txt', 'allr3n08.txt', 'allr3n09.txt',
    None, None, None]

rank_4 = [None, None, None, None, 'allr4n04.txt',
    'allr4n05.txt', 'allr4n06.txt', 'allr4n07.txt', 'allr4n08.txt', None,
    None, None, None]

extensions = [None, rank_1, rank_2, rank_3, rank_4]

groundset = 'abcdefghijklm'

# Get a list of non-isomorphic matroids of rank r on
# n elements from Hiraishi's Database
def get_matroids_database(r,n):
    matroids = []
    if (r <= 4 and n <= 12):
        if not extensions[r][n]:
            return matroids

        # Open the relavent page in the database
        link = url_base + extensions[r][n]
        f = urlopen(link)
        myfile = f.read().decode('utf-8')
        # Split the html into RevLex lines
        lines = myfile.split('\r\n')

        # convert each nonzero RevLex into a matroid
        for line in lines:
            if line == '':
                continue

            M = Matroid(groundset[:n], line, rank=r)
            if M.is_valid():
                matroids.append(M)

    return matroids

rank_1_simple = [None, 'simpler1n01.txt', None, None, None, None, None, None, None,
    None, None, None, None]

rank_2_simple = [None, None, 'simpler2n02.txt', 'simpler2n03.txt', 'simpler2n04.txt',
    'simpler2n05.txt', 'simpler2n06.txt', 'simpler2n07.txt', 'simpler2n08.txt',
    'simpler2n09.txt', 'simpler2n10.txt', 'simpler2n11.txt', 'simpler2n12.txt']

rank_3_simple = [None, None, None, 'simpler3n03.txt', 'simpler3n04.txt', 'simpler3n05.txt',
    'simpler3n06.txt', 'simpler3n07.txt', 'simpler3n08.txt', 'simpler3n09.txt', 'simpler3n10.txt',
    None, None]

rank_4_simple = [None, None, None, None, 'simpler4n04.txt', 'simpler4n05.txt', 'simpler4n06.txt',
    'simpler4n07.txt', 'simpler4n08.txt', None, None, None, None]

simple_extensions = [None, rank_1_simple, rank_2_simple, rank_3_simple, rank_4_simple]

# Get a list of non-isomorphic simple matroids of rank r on
# n elements from Hiraishi's Database
def get_simple_matroids_database(r,n):
    matroids = []
    if (r <= 4 and n <= 12):
        if not simple_extensions[r][n]:
            return matroids

        # Open the relavent page in the database
        link = url_base + simple_extensions[r][n]
        f = urlopen(link)
        myfile = f.read().decode('utf-8')
        # Split the html into RevLex lines
        lines = myfile.split('\n')

        # convert each nonzero RevLex into a matroid
        for line in lines:
            if line == '':
                continue

            M = Matroid(groundset[:n], line, rank=r)
            if M.is_valid():
                matroids.append(M)

    return matroids