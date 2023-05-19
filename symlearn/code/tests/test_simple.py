import os
import subprocess


def learn_ft(file, *args):
    file = os.path.join("../examples/simple_tests", file)
    output = subprocess.check_output(["python3", "learn_ft.py", file, '-d'] + list(args), stderr=subprocess.PIPE)
    ft = None
    eval_c = 0
    eval_d = 0
    eval_s = 0
    for line in output.decode('utf-8').split('\n'):
        if line.startswith("Learned fault tree:"):
            ft = line[len("Learned fault tree: "):].strip()
        if line.startswith("Error w.r.t. cut sets (ϕ_c):"):
            eval_c = float(line[len("Error w.r.t. cut sets (ϕ_c):"):].strip())
        if line.startswith("Error w.r.t. data set (ϕ_d):"):
            eval_d = float(line[len("Error w.r.t. data set (ϕ_d):"):].strip())
        if line.startswith("Size of fault tree (ϕ_s):"):
            eval_s = float(line[len("Size of fault tree (ϕ_s):"):].strip())

    return ft, eval_c, eval_d, eval_s


def test_simple():
    ft, eval_c, eval_d, eval_s = learn_ft("simple.csv")
    assert eval_c == 0.0
    assert eval_d == 0.0
    assert ft == "OR(A, AND(B, C))"


def test_simple_symmetry_or():
    ft, eval_c, eval_d, eval_s = learn_ft("simple_symmetry_or.csv")
    assert eval_c == 0.0
    assert eval_d == 0.0
    assert ft == "OR(AND(A, OR(B, C)), AND(D, OR(E, F)))"


def test_simple_symmetry_or_mult():
    ft, eval_c, eval_d, eval_s = learn_ft("simple_symmetry_or_mult.csv")
    assert eval_c == 0.0
    assert eval_d == 0.0
    assert ft == "OR(AND(E, OR(F, H, G)), AND(A, OR(B, D, C)))"


def test_simple_symmetry_or3():
    ft, eval_c, eval_d, eval_s = learn_ft("simple_symmetry_or3.csv")
    assert eval_c == 0.0
    assert eval_d == 0.0
    assert ft == "OR(AND(A, OR(B, C)), AND(G, OR(H, I)), AND(D, OR(E, F)))"


def test_simple_symmetry_and():
    ft, eval_c, eval_d, eval_s = learn_ft("simple_symmetry_and.csv")
    assert eval_c == 0.0
    assert eval_d == 0.0
    assert ft == "AND(AND(A, D), OR(AND(E, OR(B, C)), AND(F, OR(C, B))))"


def test_simple_symmetry_shared_or():
    ft, eval_c, eval_d, eval_s = learn_ft("simple_symmetry_shared_or.csv")
    assert eval_c == 0.0
    assert eval_d == 0.0
    assert ft == "OR(AND(C, OR(A, B)), AND(E, OR(D, B)))"


def test_simple_symmetry_shared_and():
    ft, eval_c, eval_d, eval_s = learn_ft("simple_symmetry_shared_and.csv")
    assert eval_c == 0.0
    assert eval_d == 0.0
    assert ft == "AND(OR(C, AND(B, A)), OR(E, AND(B, D)))"
