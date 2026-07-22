from scripts.generate_api_reference import extract_functions_from_file


def test_generator_indexes_importable_functions_with_complete_signature(tmp_path):
    module = tmp_path / "example.py"
    module.write_text(
        """
def public(a: int, /, b=1, *values: float, required, option=None, **kwargs) -> str:
    \"\"\"Public summary.\"\"\"
    return \"ok\"

class Example:
    def method(self):
        return None

def outer():
    def nested():
        return None
    return nested()
""",
        encoding="utf-8",
    )

    functions = extract_functions_from_file(module, "example", "hic_basic/example.py")

    assert [function["name"] for function in functions] == ["public", "outer"]
    assert functions[0]["signature"] == (
        "public(a: int, /, b=1, *values: float, required, option=None, **kwargs)"
    )
