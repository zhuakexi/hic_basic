"use strict";

const fs = require("fs");
const path = require("path");
const { pathToFileURL } = require("url");

async function main() {
  if (process.argv.length !== 3) throw new Error("render driver requires one request JSON path");
  const request = JSON.parse(fs.readFileSync(process.argv[2], "utf8"));
  const playwrightModule = process.env.HIC_BASIC_PLAYWRIGHT_MODULE;
  if (!playwrightModule) throw new Error("HIC_BASIC_PLAYWRIGHT_MODULE is not set");
  const { chromium } = require(playwrightModule);
  const launchOptions = {
    headless: true,
    args: request.browser.args || [],
  };
  if (request.browser.executable_path) launchOptions.executablePath = request.browser.executable_path;

  const browser = await chromium.launch(launchOptions);
  try {
    const context = await browser.newContext({
      acceptDownloads: true,
      viewport: {
        width: request.browser.viewport[0],
        height: request.browser.viewport[1],
      },
    });
    context.setDefaultTimeout(request.browser.timeout_ms);
    await context.route("**/*", async (route) => {
      const url = route.request().url();
      if (url.startsWith("http://") || url.startsWith("https://")) await route.abort("blockedbyclient");
      else await route.continue();
    });

    const page = await context.newPage();
    await page.goto(pathToFileURL(path.resolve(request.html_path)).href, { waitUntil: "load" });
    await page.waitForFunction(() => document.body.dataset.ready === "1");
    const protocol = await page.evaluate(() => {
      const api = window.__hicBasicRenderer;
      return api ? {
        protocolVersion: api.protocolVersion,
        renderConfigSchemaVersion: api.renderConfigSchemaVersion,
      } : null;
    });
    if (!protocol || protocol.protocolVersion !== request.protocol_version) {
      throw new Error(
        `HTML renderer protocol mismatch: expected ${request.protocol_version}, got ${protocol && protocol.protocolVersion}`,
      );
    }

    const [download, metadata] = await Promise.all([
      page.waitForEvent("download"),
      page.evaluate((config) => window.__hicBasicRenderer.downloadPng(config), request.render_config),
    ]);
    await download.saveAs(request.png_path);
    metadata.backend = {
      browserVersion: browser.version(),
      playwrightVersion: require(path.join(playwrightModule, "package.json")).version,
      protocolVersion: protocol.protocolVersion,
      renderConfigSchemaVersion: protocol.renderConfigSchemaVersion,
    };
    process.stdout.write(`${JSON.stringify(metadata)}\n`);
  } finally {
    await browser.close();
  }
}

main().catch((error) => {
  process.stderr.write(`${error.stack || error.message}\n`);
  process.exitCode = 1;
});
