#!/usr/bin/env python3

import pytest

from instanexus.preprocessing import remove_modifications

def test_remove_modifications():
    assert remove_modifications("A(ox)BC(mod)D") == "ABCD"
    assert remove_modifications("A[UNIMOD:21]BC[UNIMOD:35]D") == "ABCD"
    assert (
        remove_modifications("A(ox)[UNIMOD:21]BC(mod)[UNIMOD:35]D") == "ABCD"
    )
    assert remove_modifications(None) is None
    assert remove_modifications("ACD") == "ACD"
    assert remove_modifications("A(I)BCD") == "ABCD"
    assert remove_modifications("A(ox)B(I)C(mod)D") == "ABCD"
    assert (
        remove_modifications("A(ox)[UNIMOD:21]B(I)C(mod)[UNIMOD:35]D") == "ABCD"
    )
    assert remove_modifications("AI BCD") == "AL BCD"
    assert remove_modifications("A(ox)I B(mod)CD") == "AL BCD"