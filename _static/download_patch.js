// Adapted from TeachBooks/JupyterBook-Patches:
// https://github.com/TeachBooks/JupyterBook-Patches/blob/main/src/jupyterbook_patches/patches/_static/download_patch.js
// Download notebook and Markdown sources instead of opening them in the browser.
document.addEventListener('DOMContentLoaded', function () {
    const dropdown = document.querySelector('.dropdown-download-buttons');
    if (dropdown) {
        dropdown.querySelectorAll('a').forEach(a => {
            const href = a.getAttribute('href');
            if (href && (href.endsWith('.md') || href.endsWith('.ipynb'))) {
                if (!a.hasAttribute('download')) {
                    a.setAttribute('download', '');
                }
            }
        });
    }
});
