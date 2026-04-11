#!/usr/bin/env python3
"""Add copy-to-clipboard button to all step HTML files.

Inserts:
1. CSS for .code-block-wrapper and .copy-btn into <style>
2. JS at end of <body> to wrap code blocks and add copy buttons
"""
import glob, re, os

STEPS_DIR = os.path.join(os.path.dirname(__file__), '..', 'steps')

# CSS to add (right after existing .code-block rule)
COPY_CSS = """        .code-block-wrapper {position: relative;}
        .copy-btn {position: absolute; top: 6px; right: 6px; background: rgba(255,255,255,0.15); color: rgba(255,255,255,0.7); border: 1px solid rgba(255,255,255,0.2); border-radius: 4px; padding: 3px 10px; font-size: 0.75em; cursor: pointer; font-family: 'Segoe UI', sans-serif; transition: all 0.2s; z-index: 2; line-height: 1.4;}
        .copy-btn:hover {background: rgba(255,255,255,0.25); color: white;}
        .copy-btn.copied {background: rgba(76,175,80,0.4); color: #c8e6c9; border-color: rgba(76,175,80,0.5);}"""

# JS to add before closing </body>
COPY_JS = """    <!-- Copy-to-clipboard for code blocks -->
    <script>
    (function() {
        document.querySelectorAll('.code-block').forEach(function(block) {
            /* Skip if already inside wrapper or is a data display (scrollable sample tables) */
            if (block.parentElement.classList.contains('code-block-wrapper')) return;
            var wrapper = document.createElement('div');
            wrapper.className = 'code-block-wrapper';
            block.parentNode.insertBefore(wrapper, block);
            wrapper.appendChild(block);
            var btn = document.createElement('button');
            btn.className = 'copy-btn';
            btn.textContent = 'Copy';
            btn.setAttribute('aria-label', 'Copy code to clipboard');
            btn.addEventListener('click', function() {
                var text = block.textContent;
                /* Strip leading $ prompts from command lines */
                text = text.replace(/^\\$ /gm, '');
                navigator.clipboard.writeText(text).then(function() {
                    btn.textContent = 'Copied!';
                    btn.classList.add('copied');
                    setTimeout(function() {
                        btn.textContent = 'Copy';
                        btn.classList.remove('copied');
                    }, 2000);
                });
            });
            wrapper.appendChild(btn);
        });
    })();
    </script>"""


def patch_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    modified = False

    # 1. Check if already patched
    if '.copy-btn' in content:
        return False

    # 2. Add CSS after .code-block rule
    # Find the .code-block CSS line
    css_pattern = r'(\.code-block \{[^}]+\})'
    m = re.search(css_pattern, content)
    if m:
        insert_pos = m.end()
        content = content[:insert_pos] + '\n' + COPY_CSS + content[insert_pos:]
        modified = True
    else:
        print(f"  WARNING: no .code-block CSS found in {os.path.basename(filepath)}")

    # 3. Add JS before </body>
    body_close = content.rfind('</body>')
    if body_close >= 0:
        content = content[:body_close] + COPY_JS + '\n' + content[body_close:]
        modified = True
    else:
        print(f"  WARNING: no </body> found in {os.path.basename(filepath)}")

    if modified:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
    return modified


def main():
    files = sorted(glob.glob(os.path.join(STEPS_DIR, 'step*.html')))
    files += glob.glob(os.path.join(STEPS_DIR, 'next_steps.html'))
    # Exclude _old files
    files = [f for f in files if '_old' not in os.path.basename(f)]

    print(f"Found {len(files)} files to patch")
    patched = 0
    for f in sorted(files):
        name = os.path.basename(f)
        result = patch_file(f)
        status = "PATCHED" if result else "SKIPPED (already has .copy-btn)"
        print(f"  {name}: {status}")
        if result:
            patched += 1
    print(f"\nDone: {patched} files patched, {len(files) - patched} skipped")


if __name__ == '__main__':
    main()
