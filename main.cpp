#include <fstream>
#include <iostream>

#include <stdexcept>

#include <cmath>

#include <map>
#include <set>
#include <string>
#include <vector>

using namespace std;

const size_t bloomSize = 480000;
const size_t bloomK = 1;

const string bigramFrequency[] = { "jq", "jx", "jz", "qj", "qx", "qz", "vq", "xj", "zx", "gx", "hx", "vj", "wx", "qg", "qh", "qp", "qv", "qy", "wq", "kx", "qk", "vb", "zj", "fq", "px", "qc", "qf", "vw", "vx", "xz", "cx", "fx", "jv", "kz", "qn", "xk", "cj", "gq", "jb", "kq", "qb", "bx", "jg", "pq", "qd", "qm", "qw", "sx", "fv", "fz", "jf", "jw", "mx", "qq", "tx", "dx", "ql", "qt", "vf", "pz", "vz", "jm", "qo", "wj", "zq", "bq", "jt", "qr", "yq", "jj", "qe", "zf", "jl", "mq", "vh", "jh", "vp", "jp", "xq", "bz", "xg", "fj", "xx", "jk", "jy", "xr", "vg", "tq", "dq", "jc", "qs", "vm", "lx", "xn", "yj", "zg", "jd", "xv", "rx", "wv", "cv", "hq", "jr", "xd", "mz", "yy", "gz", "js", "vc", "zp", "gv", "hj", "jn", "zv", "hz", "lq", "gj", "mj", "pv", "wz", "vk", "vt", "pj", "vd", "cw", "zr", "zd", "zc", "fk", "fg", "zt", "qi", "xw", "zn", "zw", "fw", "lj", "zs", "cf", "fp", "kj", "vn", "zk", "fh", "xm", "bk", "kg", "xf", "fd", "xb", "cb", "mk", "fb", "zm", "wg", "cg", "qa", "zb", "vv", "gk", "dk", "cp", "fc", "uw", "xl", "mg", "uq", "gc", "sj", "fm", "bg", "pg", "sz", "cz", "pk", "ww", "vl", "zh", "fn", "uu", "tj", "dz", "hv", "gp", "lz", "kc", "kd", "kv", "wc", "mv", "xs", "bw", "hg", "nx", "mh", "tk", "pd", "cm", "cd", "uj", "hh", "mw", "vs", "wp", "bv", "hk", "bf", "mt", "gf", "yz", "yk", "md", "wf", "wu", "vr", "kp", "tv", "hd", "bp", "hc", "rq", "xh", "bn", "vy", "uh", "gd", "pc", "iy", "pw", "wm", "rj", "iw", "yv", "pf", "hp", "wk", "oj", "cq", "ij", "kk", "mr", "wt", "yx", "yh", "pb", "td", "kf", "rz", "uy", "dp", "bh", "yf", "zu", "pm", "bj", "gb", "bm", "hf", "lr", "dj", "wy", "mf", "dc", "xu", "dt", "zl", "wb", "wd", "gw", "cn", "kb", "km", "tg", "sg", "oq", "kt", "sr", "lh", "yu", "yw", "kw", "gt", "bc", "hb", "df", "sd", "uz", "ej", "lw", "fs", "nq", "aq", "ml", "bd", "bt", "dv", "tp", "vu", "zy", "ih", "aj", "ux", "mc", "pn", "db", "uv", "xc", "aa", "hw", "ao", "iq", "nj", "sf", "ji", "tb", "dh", "dw", "yb", "kr", "ez", "nz", "ii", "wl", "xy", "tf", "tz", "fy", "zz", "xo", "dm", "sv", "wr", "kh", "by", "ku", "ky", "gm", "uk", "xp", "ws", "ix", "eq", "kn", "tn", "lg", "oz", "yg", "nw", "uo", "tm", "ek", "mn", "sq", "rw", "xa", "sb", "lk", "xt", "uf", "lf", "lb", "hm", "xe", "dg", "oh", "cs", "hs", "tw", "lv", "hn", "ln", "oy", "ft", "nm", "nh", "ax", "je", "wn", "jo", "yi", "lc", "kl", "lp", "nb", "yd", "ju", "nl", "wh", "rf", "bs", "dn", "eh", "sw", "nr", "py", "lm", "nv", "ko", "az", "zi", "yo", "bb", "np", "dy", "rv", "hl", "ny", "gy", "ox", "za", "zo", "tc", "ye", "ja", "xi", "rh", "dd", "yr", "gg", "of", "ok", "ks", "yt", "aw", "rk", "af", "yc", "ym", "ah", "cc", "sk", "nk", "vo", "ug", "ik", "ht", "gn", "ya", "dl", "my", "iu", "ew", "yn", "sy", "ey", "nf", "sn", "oe", "ze", "gs", "rl", "fr", "iz", "tl", "rb", "ld", "gh", "ff", "ms", "nu", "rp", "sl", "eb", "ds", "cy", "ay", "yp", "hu", "fu", "ef", "mm", "pp", "du", "ak", "ud", "ka", "ls", "ys", "ev", "ue", "oa", "lt", "pu", "av", "wo", "yl", "fa", "ps", "fl", "hr", "rg", "wi", "we", "gu", "ua", "ui", "pt", "mu", "eg", "mb", "ib", "eu", "ex", "uc", "up", "rn", "ub", "ob", "eo", "gl", "fe", "nn", "ei", "ow", "cl", "dr", "go", "fo", "if", "oi", "va", "bu", "rc", "ki", "qu", "wa", "ov", "rr", "rm", "au", "rd", "mp", "ru", "br", "lu", "pl", "fi", "ip", "ty", "tt", "cu", "ig", "iv", "tu", "ry", "hy", "ai", "ts", "ae", "gr", "ck", "bo", "ir", "od", "ep", "ee", "cr", "oo", "sp", "sm", "im", "gi", "ke", "ga", "ag", "vi", "ut", "um", "do", "so", "rt", "ct", "bi", "su", "sc", "be", "ba", "og", "pi", "ad", "bl", "ap", "da", "oc", "ci", "sa", "em", "ab", "sh", "ec", "ge", "po", "mo", "pa", "rs", "nc", "am", "pr", "ot", "ul", "op", "id", "ce", "ly", "ur", "nd", "ea", "os", "ho", "hi", "ns", "ve", "ph", "ie", "om", "ha", "pe", "th", "ou", "mi", "ac", "si", "no", "et", "di", "ol", "as", "il", "us", "lo", "na", "me", "tr", "he", "ll", "io", "ia", "ta", "el", "ch", "ma", "ca", "to", "de", "un", "ni", "ss", "co", "it", "ed", "se", "la", "nt", "ng", "or", "ro", "ic", "li", "ar", "st", "ne", "ra", "ri", "le", "al", "re", "at", "en", "is", "te", "an", "ti", "on", "es", "in", "er"};
const size_t bigramCount = sizeof(bigramFrequency) / sizeof(bigramFrequency[0]);
const size_t strongBigramUsed = 40;
const size_t bigramUsed = 80;

const set<string> disallowedBegins = {
"bq", "fk", "gx", "hx", "jq", "jx", "jz", "kx", "kz", "lk", "lq", "mq", "qg", "qj", "qx", "qz", "rz", "uo", "uq", "vk", "vq", "vz", "wq", "wx", "wz", "xg", "xj", "xk", "xz", "yj", "yk", "yx", "yz", "zc", "zf", "zj", "zq", "zv", "zx" };

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
size_t countConsolantSeq(const string& word) {
    size_t bestCount = 0;
    size_t n = word.size();
    for (size_t i = 0; i < n; ++i) {
        size_t j = i;
        while (j < n && testVowel(word[j]))
            ++j;
        bestCount = max(bestCount, j - i);
    }
    return bestCount;
}
size_t countVowelSeq(const string& word) {
    size_t bestCount = 0;
    size_t n = word.size();
    for (size_t i = 0; i < n; ++i) {
        size_t j = i;
        while (j < n && !testVowel(word[j]))
            ++j;
        bestCount = max(bestCount, j - i);
    }
    return bestCount;
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

map<string, size_t> initBigramMap() {
    map<string, size_t> bigramMap;
    for (size_t i = 0; i < bigramCount; ++i) {
        bigramMap[bigramFrequency[i]] = i;
    }
    return bigramMap;
}

bool testWord(const string& word, const vector<bool>& bloom) {
    size_t n = word.size();
    string wordForBloom = word;
    bool isSword = false;

    if (n > 1) {
        string prefix = word.substr(0, 2);
        if (disallowedBegins.find(prefix) != disallowedBegins.end())
            return false;
    }
    if (n > 2 && word[n - 1] == 's' && word[n-2] == '\'') {
        wordForBloom = word.substr(0, n - 2);
        isSword = true;
    }
    if (!testBloom(wordForBloom, bloom))
        return false;

    if (word != wordForBloom)
        if (!testWord(wordForBloom, bloom))
            return false;

    size_t m = isSword
        ? n - 2
        : n;
    size_t containAp = countAp(word);
    size_t strangeCount = 0;
    if (m < 3) {
        for (size_t k = 0; k < strongBigramUsed; ++k) {
            if (!testBigram(word, bigramFrequency[k]))
                return false;
        }
    } else {
        for (size_t k = 0; k < bigramUsed; ++k) {
            if (!testBigram(word, bigramFrequency[k]))
                return false;
        }
    }
    static map<string, size_t> bigramMap = initBigramMap();
    size_t bigramPoints = 0;
    for (size_t i = 1; i < m; ++i) {
        size_t bigramPoint = bigramMap[word.substr(i-1, 2)];
        bigramPoints += bigramPoint * bigramPoint;
    }
    double bigramPointCoeff[] = {
        0.0,
        0.0,
        0.0,
        0.24,
        0.24,
        0.24,
        0.24,
        0.24,
        0.24,
        0.24,
        0.23, // m == 10
        0.23,
        0.23,
        0.23,
        0.23,
        0.23,
        0.23,
        0.23,
        0.23,
        0.23,
        0.23, // m == 20
        0.23,
        0.23,
    };
    if (sqrt(bigramPoints) < bigramPointCoeff[m] * bigramCount * (m - 1)) {
        return false;
    }
    if (containAp && !isSword)
        return false;
    if (n < 3 && containAp)
        return false;
    if (m < 3)
        return true;

    if (testTripple(word) && n != 3)
        return false;

    if (strangeCount > 1)
        return false;

    if (n > 16) {
        // too many correct words is lost, but trash filtered better (at this moment).
        return false;
    }

    size_t consolantSeq = countConsolantSeq(word);
    size_t vowelSeq = countVowelSeq(word);
    size_t consolantCount = countConsolant(word);
    size_t vowelCount = n - consolantCount;
    if (consolantSeq > 6)
        return false;
    if (vowelSeq > 6)
        return false;
    if (vowelCount == 1 && endedByEd(word))
        ++strangeCount;
    if (m > 3 && vowelCount == 0)
        return false;
    if (m > 3 && consolantCount == 0)
        return false;
    if (n > 7 && vowelCount < n / 5)
        ++strangeCount;
//    if (n > 7 && consolantCount < n / 4)
 //       ++strangeCount;
    if (strangeCount > 0)
        return false;

//    if (word[0] >= 'z')
 //       ++strangeCount;

  //  if (countPair(word) > 3)
   //     ++strangeCount;
 //   if (strangeCount > 1)
  //      return false;

//    if (containAp)
//        return false;
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
        cerr << "usage: " << argv[0] << " DIR" << endl;
        return 1;
    }

    size_t countTrashOk = 0;
    size_t countTrashFail = 0;
    size_t countWordsOk = 0;
    size_t countWordsFail = 0;

    string lex;

    string dictionaryFilename = "words-2.txt";
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
            if (rand() < 800000)
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
            if (rand() < 2000000)
                cerr << lex << endl;
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
    return 0;

    const set<string> allowed = {
"aa", "ab", "ac", "ad", "ae", "af", "ag", "ah", "ai", "aj", "ak", "al", "am", "an", "ao", "ap", "aq", "ar", "as", "at", "au", "av", "aw", "ax", "ay", "az", "ba", "bb", "bc", "bd", "be", "bf", "bg", "bh", "bi", "bj", "bk", "bl", "bm", "bn", "bo", "bp", "br", "bs", "bt", "bu", "bv", "bw", "bx", "by", "bz", "ca", "cb", "cc", "cd", "ce", "cf", "cg", "ch", "ci", "cj", "ck", "cl", "cm", "cn", "co", "cp", "cq", "cr", "cs", "ct", "cu", "cv", "cw", "cx", "cy", "cz", "da", "db", "dc", "dd", "de", "df", "dg", "dh", "di", "dj", "dk", "dl", "dm", "dn", "do", "dp", "dq", "dr", "ds", "dt", "du", "dv", "dw", "dx", "dy", "dz", "ea", "eb", "ec", "ed", "ee", "ef", "eg", "eh", "ei", "ej", "ek", "el", "em", "en", "eo", "ep", "eq", "er", "es", "et", "eu", "ev", "ew", "ex", "ey", "ez", "fa", "fb", "fc", "fd", "fe", "ff", "fg", "fh", "fi", "fj", "fl", "fm", "fn", "fo", "fp", "fq", "fr", "fs", "ft", "fu", "fv", "fw", "fx", "fy", "fz", "ga", "gb", "gc", "gd", "ge", "gf", "gg", "gh", "gi", "gj", "gk", "gl", "gm", "gn", "go", "gp", "gq", "gr", "gs", "gt", "gu", "gv", "gw", "gy", "gz", "ha", "hb", "hc", "hd", "he", "hf", "hg", "hh", "hi", "hj", "hk", "hl", "hm", "hn", "ho", "hp", "hq", "hr", "hs", "ht", "hu", "hv", "hw", "hy", "hz", "ia", "ib", "ic", "id", "ie", "if", "ig", "ih", "ii", "ij", "ik", "il", "im", "in", "io", "ip", "iq", "ir", "is", "it", "iu", "iv", "iw", "ix", "iy", "iz", "ja", "jb", "jc", "jd", "je", "jf", "jg", "jh", "ji", "jj", "jk", "jl", "jm", "jn", "jo", "jp", "jr", "js", "jt", "ju", "jv", "jw", "jy", "ka", "kb", "kc", "kd", "ke", "kf", "kg", "kh", "ki", "kj", "kk", "kl", "km", "kn", "ko", "kp", "kq", "kr", "ks", "kt", "ku", "kv", "kw", "ky", "la", "lb", "lc", "ld", "le", "lf", "lg", "lh", "li", "lj", "ll", "lm", "ln", "lo", "lp", "lr", "ls", "lt", "lu", "lv", "lw", "lx", "ly", "lz", "ma", "mb", "mc", "md", "me", "mf", "mg", "mh", "mi", "mj", "mk", "ml", "mm", "mn", "mo", "mp", "mr", "ms", "mt", "mu", "mv", "mw", "mx", "my", "mz", "na", "nb", "nc", "nd", "ne", "nf", "ng", "nh", "ni", "nj", "nk", "nl", "nm", "nn", "no", "np", "nq", "nr", "ns", "nt", "nu", "nv", "nw", "nx", "ny", "nz", "oa", "ob", "oc", "od", "oe", "of", "og", "oh", "oi", "oj", "ok", "ol", "om", "on", "oo", "op", "oq", "or", "os", "ot", "ou", "ov", "ow", "ox", "oy", "oz", "pa", "pb", "pc", "pd", "pe", "pf", "pg", "ph", "pi", "pj", "pk", "pl", "pm", "pn", "po", "pp", "pq", "pr", "ps", "pt", "pu", "pv", "pw", "px", "py", "pz", "qa", "qb", "qc", "qd", "qe", "qf", "qh", "qi", "qk", "ql", "qm", "qn", "qo", "qp", "qq", "qr", "qs", "qt", "qu", "qv", "qw", "qy", "ra", "rb", "rc", "rd", "re", "rf", "rg", "rh", "ri", "rj", "rk", "rl", "rm", "rn", "ro", "rp", "rq", "rr", "rs", "rt", "ru", "rv", "rw", "rx", "ry", "sa", "sb", "sc", "sd", "se", "sf", "sg", "sh", "si", "sj", "sk", "sl", "sm", "sn", "so", "sp", "sq", "sr", "ss", "st", "su", "sv", "sw", "sx", "sy", "sz", "ta", "tb", "tc", "td", "te", "tf", "tg", "th", "ti", "tj", "tk", "tl", "tm", "tn", "to", "tp", "tq", "tr", "ts", "tt", "tu", "tv", "tw", "tx", "ty", "tz", "ua", "ub", "uc", "ud", "ue", "uf", "ug", "uh", "ui", "uj", "uk", "ul", "um", "un", "up", "ur", "us", "ut", "uu", "uv", "uw", "ux", "uy", "uz", "va", "vb", "vc", "vd", "ve", "vf", "vg", "vh", "vi", "vj", "vl", "vm", "vn", "vo", "vp", "vr", "vs", "vt", "vu", "vv", "vw", "vx", "vy", "wa", "wb", "wc", "wd", "we", "wf", "wg", "wh", "wi", "wj", "wk", "wl", "wm", "wn", "wo", "wp", "wr", "ws", "wt", "wu", "wv", "ww", "wy", "xa", "xb", "xc", "xd", "xe", "xf", "xh", "xi", "xl", "xm", "xn", "xo", "xp", "xq", "xr", "xs", "xt", "xu", "xv", "xw", "xx", "xy", "ya", "yb", "yc", "yd", "ye", "yf", "yg", "yh", "yi", "yl", "ym", "yn", "yo", "yp", "yq", "yr", "ys", "yt", "yu", "yv", "yw", "yy", "za", "zb", "zd", "ze", "zg", "zh", "zi", "zk", "zl", "zm", "zn", "zo", "zp", "zr", "zs", "zt", "zu", "zw", "zy", "zz"};
    for (size_t c1 = 'a'; c1 <= 'z'; ++c1) {
        for (size_t c2 = 'a'; c2 <= 'z'; ++c2) {
            string s;
            s += c1;
            s += c2;
            if (allowed.find(s) == allowed.end())
                cout << s << endl;
        }
    }
}
