"""Convert the inline/display LaTeX used in _v3_appendix.md to Unicode text for the Word build."""
import re

CMD = {r"\alpha": "α", r"\beta": "β", r"\lambda": "λ", r"\chi": "χ", r"\nu": "ν", r"\sigma": "σ", r"\varepsilon": "ε",
       r"\mu": "μ", r"\le": "≤", r"\ge": "≥", r"\in": "∈", r"\infty": "∞", r"\pm": "±", r"\sum": "Σ", r"\dots": "…",
       r"\mid": " | ", r"\lceil": "⌈", r"\rceil": "⌉", r"\{": "{", r"\}": "}", r"\#": "#", r"\Big": "", r"\ ": " ",
       r"\quad": " ", r"\max": "max", r"{:}": ":"}
SUB = dict(zip("0123456789+-=()aehijklmnoprstuvxβν", "₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ₐₑₕᵢⱼₖₗₘₙₒₚᵣₛₜᵤᵥₓᵦᵥ"))
SUP = dict(zip("0123456789+-=()n", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾ⁿ"))

def script(s, table, mark):
    return ("\x01" if mark == "_" else "\x03") + s + ("\x02" if mark == "_" else "\x04")

def conv(t):
    t = re.sub(r"\\text\{([^}]*)\}", r"\1", t)
    t = re.sub(r"\{\\rm ([^}]*)\}", r"\1", t)
    t = t.replace(r"\hat\mu", "μ\u0302").replace(r"\hat q", "q\u0302")
    for k in sorted(CMD, key=len, reverse=True):
        t = t.replace(k, CMD[k])
    for _ in range(3):
        t = re.sub(r"_\{([^{}]*)\}", lambda m: script(m.group(1).replace(" ", ""), SUB, "_"), t)
        t = re.sub(r"\^\{([^{}]*)\}", lambda m: script(m.group(1).replace(" ", ""), SUP, "^"), t)
    t = re.sub(r"_(\w)", lambda m: script(m.group(1), SUB, "_"), t)
    t = re.sub(r"\^(\w)", lambda m: script(m.group(1), SUP, "^"), t)
    t = t.replace("\x01", "_{").replace("\x02", "}").replace("\x03", "^{").replace("\x04", "}")
    return t.replace("-", "−").replace("  ", " ")

def convert(md):
    md = re.sub(r"\$\$\s*(.*?)\s*\$\$", lambda m: "\n<center>" + conv(m.group(1)) + "\n", md, flags=re.S)
    return re.sub(r"\$([^$]+)\$", lambda m: conv(m.group(1)), md)

if __name__ == "__main__":
    out = convert(open("_v3_appendix.md").read())
    print(out[:3000])
    print(sorted(set(re.findall(r"\\\w+|[_^]", out))))
