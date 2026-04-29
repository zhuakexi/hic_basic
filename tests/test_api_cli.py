import json
import subprocess
import sys


def run_api(*args):
    return subprocess.run(
        [sys.executable, "-m", "hic_basic", "api", *args],
        capture_output=True,
        text=True,
        check=False,
    )


def test_api_list_module_contains_cool2mat():
    result = run_api("list", "coolstuff")

    assert result.returncode == 0, result.stderr
    assert "cool2mat" in result.stdout


def test_api_search_json_is_parseable():
    result = run_api("search", "cool2mat", "--json")

    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["query"] == "cool2mat"
    assert payload["matches"]
    assert payload["matches"][0]["symbol"].startswith("hic_basic.")


def test_api_show_json_has_core_fields():
    result = run_api("show", "cool2mat", "--json")

    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["symbol"] == "hic_basic.coolstuff.cool2mat"
    assert payload["signature"].startswith("cool2mat(")
    assert payload["import"] == "from hic_basic.coolstuff import cool2mat"


def test_api_where_reports_source():
    result = run_api("where", "cool2mat")

    assert result.returncode == 0, result.stderr
    assert "hic_basic.coolstuff.cool2mat" in result.stdout
    assert ".py:" in result.stdout


def test_api_show_ambiguous_symbol_lists_candidates():
    result = run_api("show", "cool")

    assert result.returncode == 0, result.stderr
    assert "Ambiguous symbol: cool" in result.stdout
    assert "hic_basic api show" in result.stdout
