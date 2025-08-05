"use client"

import Image from "next/image"
import { useState, useEffect } from "react"

export default function Header() {
  const [isVisible, setIsVisible] = useState(true)
  const [lastScrollY, setLastScrollY] = useState(0)

  useEffect(() => {
    const controlHeader = () => {
      const currentScrollY = window.scrollY

      // Show header when scrolling up, hide when scrolling down
      if (currentScrollY > lastScrollY && currentScrollY > 100) {
        setIsVisible(false)
      } else {
        setIsVisible(true)
      }

      setLastScrollY(currentScrollY)
    }

    window.addEventListener("scroll", controlHeader)

    return () => {
      window.removeEventListener("scroll", controlHeader)
    }
  }, [lastScrollY])

  return (
    <header
      className={`fixed top-0 left-0 w-full bg-background text-foreground shadow-md z-50 flex items-center justify-between px-4 py-2 transition-transform duration-300 ease-in-out ${
        isVisible ? "translate-y-0" : "-translate-y-full"
      }`}
    >
      <div className="flex items-center" style={{ paddingLeft: "2rem" }}>
        <Image src="/aganitha-logo.png" alt="Aganitha Logo" width={140} height={150} style={{ objectFit: "contain" }} />
      </div>
      {/* Centered Content */}
      <div className="flex flex-col items-center">
        <h1 className="text-2xl font-bold text-primary">Nitrosamine Yield Prediction</h1>
        {/* <p className="text-l text-muted-foreground">
          Discover molecular analogs by intelligently exploring the chemical neighborhoods around input compounds
        </p> */}
      </div>
      {/* Empty div to maintain balance */}
      <div className="w-12"></div> {/* Optional, to keep other elements in place if needed */}
    </header>
  )
}
