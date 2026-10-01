"use strict";

// The documentation pages: the theme switch, shared with the calculator
// (localStorage "ds-theme"), and a contents list for long pages.
(function () {
  const root = document.documentElement;

  function applyTheme(theme) {
    if (theme === "light") root.setAttribute("data-theme", "light");
    else root.removeAttribute("data-theme");
  }

  const toggle = document.getElementById("themeToggle");
  if (toggle) {
    toggle.addEventListener("click", () => {
      const theme = root.getAttribute("data-theme") === "light" ? "dark" : "light";
      applyTheme(theme);
      try { localStorage.setItem("ds-theme", theme); } catch (_) {}
    });
  }
  // A theme switched in another tab, such as the calculator, applies here too.
  window.addEventListener("storage", (event) => {
    if (event.key === "ds-theme") applyTheme(event.newValue);
  });

  // Contents: the sections of a page with four or more, shown beside the
  // article on wide screens, with the section in view marked.
  const page = document.querySelector(".page");
  const main = document.querySelector(".doc main");
  if (!page || !main) return;
  const sections = Array.from(main.children).filter((el) => el.tagName === "SECTION" && el.querySelector("h2"));
  if (sections.length < 4) return;

  const toc = document.createElement("nav");
  toc.className = "toc";
  toc.setAttribute("aria-label", "On this page");
  const title = document.createElement("p");
  title.className = "toc-title";
  title.textContent = "On this page";
  const list = document.createElement("ol");
  toc.append(title, list);

  const links = new Map();
  for (const section of sections) {
    const text = section.querySelector("h2").textContent.trim();
    if (!section.id) {
      const base = text.toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-+|-+$/g, "") || "section";
      let id = base;
      for (let n = 2; document.getElementById(id); n++) id = `${base}-${n}`;
      section.id = id;
    }
    const item = document.createElement("li");
    const link = document.createElement("a");
    link.href = `#${section.id}`;
    link.textContent = text;
    item.append(link);
    list.append(item);
    links.set(section, link);
  }
  page.append(toc);
  page.classList.add("has-toc");

  if (!("IntersectionObserver" in window)) return;
  const visible = new Set();
  const observer = new IntersectionObserver((entries) => {
    for (const entry of entries) {
      if (entry.isIntersecting) visible.add(entry.target);
      else visible.delete(entry.target);
    }
    const current = sections.find((section) => visible.has(section));
    if (!current) return;
    for (const [section, link] of links) link.classList.toggle("on", section === current);
  }, { rootMargin: "-72px 0px -60% 0px" });
  for (const section of sections) observer.observe(section);
})();
