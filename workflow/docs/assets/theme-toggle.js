(function () {
  var storageKey = "workflow-docs-theme";

  function preferredTheme() {
    var saved = localStorage.getItem(storageKey);
    if (saved === "dark" || saved === "light") return saved;
    return window.matchMedia &&
      window.matchMedia("(prefers-color-scheme: dark)").matches
      ? "dark"
      : "light";
  }

  function applyTheme(theme) {
    document.body.classList.toggle("dark-mode", theme === "dark");
    var buttons = document.querySelectorAll(".theme-toggle");
    buttons.forEach(function (button) {
      button.textContent = theme === "dark" ? "Light" : "Dark";
      button.setAttribute(
        "aria-label",
        theme === "dark" ? "Switch to light mode" : "Switch to dark mode"
      );
      button.setAttribute("title", button.getAttribute("aria-label"));
    });
  }

  function toggleTheme() {
    var next = document.body.classList.contains("dark-mode") ? "light" : "dark";
    localStorage.setItem(storageKey, next);
    applyTheme(next);
  }

  function addToggle(container) {
    if (!container || container.querySelector(".theme-toggle")) return;
    var button = document.createElement("button");
    button.className = "theme-toggle";
    button.type = "button";
    button.addEventListener("click", toggleTheme);
    container.appendChild(button);
  }

  function init() {
    applyTheme(preferredTheme());
    addToggle(document.querySelector(".wy-nav-top"));
    addToggle(document.querySelector(".wy-side-nav-search"));
    applyTheme(preferredTheme());
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
