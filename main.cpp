#include <fstream>
#include <iostream>

#include <stdexcept>

#include <map>
#include <set>
#include <string>
#include <vector>

using namespace std;

const size_t bloomSize = 480000;
const size_t bloomK = 1;

const string bigramFrequency[] = { "jq", "jx", "jz", "qj", "qx", "qz", "vq", "xj", "zx", "gx", "hx", "vj", "wx", "qg", "qh", "qp", "qv", "qy", "wq", "kx", "qk", "vb", "zj", "fq", "px", "qc", "qf", "vw", "vx", "xz", "cx", "fx", "jv", "kz", "qn", "xk", "cj", "gq", "jb", "kq", "qb", "bx", "jg", "pq", "qd", "qm", "qw", "sx", "fv", "fz", "jf", "jw", "mx", "qq", "tx", "dx", "ql", "qt", "vf", "pz", "vz", "jm", "qo", "wj", "zq", "bq", "jt", "qr", "yq", "jj", "qe", "zf", "jl", "mq", "vh", "jh", "vp", "jp", "xq", "bz", "xg", "fj", "xx", "jk", "jy", "xr", "vg", "tq", "dq", "jc", "qs", "vm", "lx", "xn", "yj", "zg", "jd", "xv", "rx", "wv", "cv", "hq", "jr", "xd", "mz", "yy", "gz", "js", "vc", "zp", "gv", "hj", "jn", "zv", "hz", "lq", "gj", "mj", "pv", "wz", "vk", "vt", "pj", "vd", "cw", "zr", "zd", "zc", "fk", "fg", "zt", "qi", "xw", "zn", "zw", "fw", "lj", "zs", "cf", "fp", "kj", "vn", "zk", "fh", "xm", "bk", "kg", "xf", "fd", "xb", "cb", "mk", "fb", "zm", "wg", "cg", "qa", "zb", "vv", "gk", "dk", "cp", "fc", "uw", "xl", "mg", "uq", "gc", "sj", "fm", "bg", "pg", "sz", "cz", "pk", "ww", "vl", "zh", "fn", "uu", "tj", "dz", "hv", "gp", "lz", "kc", "kd", "kv", "wc", "mv", "xs", "bw", "hg", "nx", "mh", "tk", "pd", "cm", "cd", "uj", "hh", "mw", "vs", "wp", "bv", "hk", "bf", "mt", "gf", "yz", "yk", "md", "wf", "wu", "vr", "kp", "tv", "hd", "bp", "hc", "rq", "xh", "bn", "vy", "uh", "gd", "pc", "iy", "pw", "wm", "rj", "iw", "yv", "pf", "hp", "wk", "oj", "cq", "ij", "kk", "mr", "wt", "yx", "yh", "pb", "td", "kf", "rz", "uy", "dp", "bh", "yf", "zu", "pm", "bj", "gb", "bm", "hf", "lr", "dj", "wy", "mf", "dc", "xu", "dt", "zl", "wb", "wd", "gw", "cn", "kb", "km", "tg", "sg", "oq", "kt", "sr", "lh", "yu", "yw", "kw", "gt", "bc", "hb", "df", "sd", "uz", "ej", "lw", "fs", "nq", "aq", "ml", "bd", "bt", "dv", "tp", "vu", "zy", "ih", "aj", "ux", "mc", "pn", "db", "uv", "xc", "aa", "hw", "ao", "iq", "nj", "sf", "ji", "tb", "dh", "dw", "yb", "kr", "ez", "nz", "ii", "wl", "xy", "tf", "tz", "fy", "zz", "xo", "dm", "sv", "wr", "kh", "by", "ku", "ky", "gm", "uk", "xp", "ws", "ix", "eq", "kn", "tn", "lg", "oz", "yg", "nw", "uo", "tm", "ek", "mn", "sq", "rw", "xa", "sb", "lk", "xt", "uf", "lf", "lb", "hm", "xe", "dg", "oh", "cs", "hs", "tw", "lv", "hn", "ln", "oy", "ft", "nm", "nh", "ax", "je", "wn", "jo", "yi", "lc", "kl", "lp", "nb", "yd", "ju", "nl", "wh", "rf", "bs", "dn", "eh", "sw", "nr", "py", "lm", "nv", "ko", "az", "zi", "yo", "bb", "np", "dy", "rv", "hl", "ny", "gy", "ox", "za", "zo", "tc", "ye", "ja", "xi", "rh", "dd", "yr", "gg", "of", "ok", "ks", "yt", "aw", "rk", "af", "yc", "ym", "ah", "cc", "sk", "nk", "vo", "ug", "ik", "ht", "gn", "ya", "dl", "my", "iu", "ew", "yn", "sy", "ey", "nf", "sn", "oe", "ze", "gs", "rl", "fr", "iz", "tl", "rb", "ld", "gh", "ff", "ms", "nu", "rp", "sl", "eb", "ds", "cy", "ay", "yp", "hu", "fu", "ef", "mm", "pp", "du", "ak", "ud", "ka", "ls", "ys", "ev", "ue", "oa", "lt", "pu", "av", "wo", "yl", "fa", "ps", "fl", "hr", "rg", "wi", "we", "gu", "ua", "ui", "pt", "mu", "eg", "mb", "ib", "eu", "ex", "uc", "up", "rn", "ub", "ob", "eo", "gl", "fe", "nn", "ei", "ow", "cl", "dr", "go", "fo", "if", "oi", "va", "bu", "rc", "ki", "qu", "wa", "ov", "rr", "rm", "au", "rd", "mp", "ru", "br", "lu", "pl", "fi", "ip", "ty", "tt", "cu", "ig", "iv", "tu", "ry", "hy", "ai", "ts", "ae", "gr", "ck", "bo", "ir", "od", "ep", "ee", "cr", "oo", "sp", "sm", "im", "gi", "ke", "ga", "ag", "vi", "ut", "um", "do", "so", "rt", "ct", "bi", "su", "sc", "be", "ba", "og", "pi", "ad", "bl", "ap", "da", "oc", "ci", "sa", "em", "ab", "sh", "ec", "ge", "po", "mo", "pa", "rs", "nc", "am", "pr", "ot", "ul", "op", "id", "ce", "ly", "ur", "nd", "ea", "os", "ho", "hi", "ns", "ve", "ph", "ie", "om", "ha", "pe", "th", "ou", "mi", "ac", "si", "no", "et", "di", "ol", "as", "il", "us", "lo", "na", "me", "tr", "he", "ll", "io", "ia", "ta", "el", "ch", "ma", "ca", "to", "de", "un", "ni", "ss", "co", "it", "ed", "se", "la", "nt", "ng", "or", "ro", "ic", "li", "ar", "st", "ne", "ra", "ri", "le", "al", "re", "at", "en", "is", "te", "an", "ti", "on", "es", "in", "er"};
const size_t bigramUsed = 80;
const size_t bigramUsedStrange = 160;

uint64_t hashString(const string& word, size_t index) {
    const uint64_t constants[] = {27644437, 115249, 33391, 108301, 115249};
    const size_t constantsCount = sizeof(constants) / sizeof(constants[0]);
    if (index + 1 > constantsCount)
        throw std::invalid_argument("overflow in hashString: index too long");

    uint64_t result = 0;
    for (auto i : word) {
        result = result * constants[index] + constants[index + 1] + i;
    }
    return result % bloomSize;
}

void setBloom(const string& word, vector<bool>& bloom) {
    for (size_t k = 0; k < bloomK; ++k) {
        bloom[hashString(word, k)] = true;
    }
}
bool testBloom(const string& word, const vector<bool>& bloom) {
    for (size_t k = 0; k < bloomK; ++k) {
        if (!bloom[hashString(word, k)])
            return false;
    }
    return true;
}

bool testVowel(char i) {
    if (i == 'a' || i == 'e' || i == 'i' || i == 'o' || i == 'u')
        return true;
    return false;
}

size_t countConsolant(const string& word) {
    size_t count = 0;
    for (auto i : word) {
        if (testVowel(i))
            ++count;
    }
    return count;
}

size_t countAp(const string& word) {
    size_t count = 0;
    for (auto i : word)
        if (i == '\'')
            ++count;
    return count;
}

bool testTripple(const string& word) {
    size_t n = word.size();
    for (size_t i = 2; i < n; ++i) {
        if (word[i] == word[i - 1] && word[i] == word[i - 2])
            return true;
    }
    return false;
}

bool testTrippleVowel(const string& word) {
    size_t n = word.size();
    for (size_t i = 2; i < n; ++i) {
        if (testVowel(word[i]) && testVowel(word[i-1]) && testVowel(word[i-2]))
            return true;
    }
    return false;
}

size_t countPair(const string& word) {
    size_t n = word.size();
    size_t count = 0;
    for (size_t i = 1; i < n; ++i) {
        if (word[i] == word[i - 1])
            ++count;
    }
    return count;
}

size_t endedByEd(const string& word) {
    size_t n = word.size();
    if (n < 4)
        return false;
    if (word[n-3] == 'd' && word[n-4] == 'e')
        return true;
    return false;
}

bool testBigram(const string& word, const string& bigram) {
    if (word.find(bigram) != string::npos)
        return false;
    return true;
}

bool testWord(const string& word, const vector<bool>& bloom) {
    if (!testBloom(word, bloom))
        return false;

    size_t n = word.size();
    size_t containAp = countAp(word);
    size_t strangeCount = 0;
    for (size_t k = 0; k < bigramUsed; ++k) {
        if (!testBigram(word, bigramFrequency[k]))
            return false;
    }
    for (size_t k = 0; k < bigramUsedStrange; ++k) {
        if (!testBigram(word, bigramFrequency[k]))
            ++strangeCount;
    }

    if (containAp > 1)
        ++strangeCount;
    if (n < 3 && containAp)
        return false;
    if (n < 3)
        return true;

    if (testTripple(word) && n != 3)
        ++strangeCount;

    if (testTrippleVowel(word) && n < 8)
        ++strangeCount;

    if (n > 15) {
        // too many correct words is lost, but trash filtered better (at this moment).
        return false;
    }
    if (strangeCount > 1)
        return false;

    if (word[n - 1] == 's' && word[n - 2] == '\'') {
        string wordBefore = word.substr(0, n - 2);
        if (!testBloom(wordBefore, bloom))
            return false;
        if (containAp == 1)
            return true;
        return false;
    }

    size_t consolantCount = countConsolant(word);
    size_t vowelCount = n - consolantCount;
    if (vowelCount == 1 && endedByEd(word))
        ++strangeCount;
    if (vowelCount == 0)
        ++strangeCount;
    if (consolantCount == 0)
        ++strangeCount;
    if (n > 7 && vowelCount < n / 5)
        ++strangeCount;
    if (n > 7 && consolantCount < n / 4)
        ++strangeCount;
    if (strangeCount > 0)
        return false;

    if (word[0] >= 'z')
        ++strangeCount;

    if (countPair(word) > 3)
        ++strangeCount;
    if (strangeCount > 1)
        return false;

    if (containAp)
        return false;
    return true;
}

void initBigram(map<string, size_t>& biMap) {
    for (char c = 'a'; c <= 'z'; ++c) {
        for (char c2 = 'a'; c2 <= 'z'; ++c2) {
            string bi;
            bi += c;
            bi += c2;
            biMap[bi] = 0;
        }
    }
}

void addBigram(const string& word, map<string, size_t>& biMap) {
    size_t n = word.size();
    for (size_t i = 1; i < n; ++i) {
        if (word[i] == '\'')
            continue;
        if (word[i-1] == '\'')
            continue;
        string bi;
        bi += word[i - 1];
        bi += word[i];
        ++biMap[bi];
    }
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        cerr << "usage: " << argv[0] << "DIR" << endl;
        return 1;
    }

    size_t countTrashOk = 0;
    size_t countTrashFail = 0;
    size_t countWordsOk = 0;
    size_t countWordsFail = 0;

    string lex;

    string dictionaryFilename = "words-uniq.txt";
    ifstream dictionaryFile(dictionaryFilename.c_str());
    set<string> dictionary;

    map<string, size_t> wordsBi;
    map<string, size_t> trashBi;
    initBigram(wordsBi);
    initBigram(trashBi);

    while (dictionaryFile >> lex) {
        addBigram(lex, wordsBi);
        dictionary.insert(lex);
    }

    vector<bool> bloom(bloomSize);
    for (auto i : dictionary) {
        setBloom(i, bloom);
    }

    string dir = argv[1];
    string trashFile = dir + "/trash.list";
    string wordsFile = dir + "/words.list";
    ifstream trash(trashFile.c_str());
    ifstream words(wordsFile.c_str());

    while (trash >> lex) {
        if (testWord(lex, bloom)) {
            ++countTrashFail;
            if (rand() < 200000)
                cerr << lex << endl;
        } else {
            ++countTrashOk;
        }
    }
    while (words >> lex) {
        addBigram(lex, wordsBi);
        if (testWord(lex, bloom)) {
            ++countWordsOk;
        } else {
            ++countWordsFail;
        }
    }

//    for (auto i : wordsBi) {
//        cout << i.second << ' ' << i.first << endl;
//    }

    double precisionTrash = 1.0 * countTrashOk / (countTrashOk + countTrashFail);
    double precisionWords = 1.0 * countWordsOk / (countWordsOk + countWordsFail);
    double precision = 1.0 * (countTrashOk + countWordsOk) / (countTrashOk + countTrashFail + countWordsFail + countWordsOk);
    cout << "Trash: " << precisionTrash << " " << countTrashOk << " : " << countTrashFail << endl;
    cout << "Words: " << precisionWords << " " << countWordsOk << " : " << countWordsFail << endl;
    cout << "Precision: " << precision << endl;
}
